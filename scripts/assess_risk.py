#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Risk Assessment Engine for AE Underreporting Evaluation.
Processes parsed lab data, CTCAE mapping, AE/MH records to classify risk levels.
"""
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


def norm(s):
    if s is None:
        return ""
    return str(s).strip()


def to_float(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    try:
        return float(norm(x))
    except (ValueError, TypeError):
        return None


def classify_visit(visit_str):
    """Classify visit phase."""
    v = norm(visit_str).upper()
    if any(k in v for k in ["筛选", "SCR", "D-", "导入期", "磨合期"]):
        return "screening"
    if re.search(r'(^|[（(])D1\b|[（(]基线[）)]|^基线$|BASELINE', v):
        return "baseline"
    return "post"


def get_visit_rank(visit_str):
    v = norm(visit_str)
    if classify_visit(v) == "screening":
        return 0
    if classify_visit(v) == "baseline":
        return 1
    m = re.search(r'D(\d+)', v)
    if m:
        return int(m.group(1))
    return 999


def is_baseline_visit(v):
    return classify_visit(v) == "baseline"


def is_screening_visit(v):
    return classify_visit(v) == "screening"


def is_abnormal(result, lo, hi):
    """Check if a result is outside normal range."""
    if result is None:
        return False, None
    if isinstance(lo, (int, float)) and result < lo:
        return True, "低"
    if isinstance(hi, (int, float)) and result > hi:
        return True, "高"
    return False, None


def get_ctcae_grade(test_name, result, lo, hi, baseline_result, ctcae_index):
    """
    Determine CTCAE grade with ULN/LLN-based calculation.
    """
    entry = ctcae_index.get(test_name, {})
    if not entry or not entry.get("has_grades"):
        return None, ""

    term = entry.get("ctcae_term", "")
    is_high = any(k in term for k in ["升高", "增加", "高", "延长", "高血压"])
    is_low = any(k in term for k in ["降低", "减少", "低", "减少"])

    if result is None:
        return None, ""

    # Weight - percentage based
    if test_name in ("体重",) and baseline_result and baseline_result != 0:
        pct = (result - baseline_result) / baseline_result * 100
        abs_pct = abs(pct)
        if abs_pct < 5:
            return None, ""
        direction = "体重增加" if pct > 0 else "体重降低"
        if abs_pct < 10:
            return 1, f"{direction}，较基线变化{pct:+.1f}%（≥5%且＜10%，1级）"
        if abs_pct < 20:
            return 2, f"{direction}，较基线变化{pct:+.1f}%（≥10%且＜20%，2级）"
        return 3, f"{direction}，较基线变化{pct:+.1f}%（≥20%，3级）"

    # QTcF
    if test_name in ("QTcF", "QT间期") and result:
        if result < 450:
            return None, ""
        if result <= 480:
            return 1, f"QTc {result}ms（450-480ms，1级）"
        if result <= 500:
            return 2, f"QTc {result}ms（481-500ms，2级）"
        excess = (result - baseline_result) if baseline_result else 0
        if result >= 501:
            level = 3
            desc = f"QTc {result}ms（≥501ms{'，较基线＞60ms' if excess > 60 else ''}，3级）"
            return level, desc
        return 2, f"QTc {result}ms"

    # For ratio-based grading (ULN/LLN)
    if is_high and hi and hi > 0:
        ratio = result / hi
        if ratio <= 1.0:
            return None, ""
        if ratio <= 3.0:
            return 1, f"＞ULN～3.0×ULN（{ratio:.1f}×ULN，1级）"
        if ratio <= 5.0:
            return 2, f"＞3.0～5.0×ULN（{ratio:.1f}×ULN，2级）"
        if ratio <= 20.0:
            return 3, f"＞5.0～20.0×ULN（{ratio:.1f}×ULN，3级）"
        return 4, f"＞20.0×ULN（{ratio:.1f}×ULN，4级）"

    if is_low and lo and lo > 0:
        ratio = result / lo
        if ratio >= 1.0:
            return None, ""
        if ratio >= 0.75:
            return 1, f"＜LLN～0.75×LLN（{ratio:.2f}×LLN，1级）"
        if ratio >= 0.5:
            return 2, f"＜0.75～0.5×LLN（{ratio:.2f}×LLN，2级）"
        if ratio >= 0.25:
            return 3, f"＜0.5～0.25×LLN（{ratio:.2f}×LLN，3级）"
        return 4, f"＜0.25×LLN（{ratio:.2f}×LLN，4级）"

    # Fallback: just check if abnormal
    if result is not None:
        if lo is not None and result < lo and is_low:
            return 1, "＜正常值下限（1级）"
        if hi is not None and result > hi and is_high:
            return 1, "＞正常值上限（1级）"

    return None, ""


def term_relevance(test_name, ae_mh_terms, direction=""):
    """
    Check if AE/MH terms can explain the lab abnormality.
    Returns (relevant: bool, reason: str).
    """
    if not ae_mh_terms:
        return False, ""

    combined_text = "；".join(ae_mh_terms)
    t = norm(test_name)

    # Urine-blood mismatch
    if any(k in combined_text for k in ("尿白细胞", "尿潜血", "尿蛋白", "尿红细胞", "尿胆原", "尿管型", "尿酮体", "尿糖")):
        if any(k in t for k in ("白细胞", "中性粒", "血红蛋白", "血小板", "ALT", "AST", "肌酐", "血糖", "胆固醇")):
            return False, "尿液检查异常与本项血液指标非同一样本系统"

    # Blood-urine mismatch
    if any(k in t for k in ("尿", "URI")):
        if any(k in combined_text for k in ("血白细胞", "贫血", "血小板", "转氨酶")):
            return False, "血液系统异常不能直接解释尿液检查异常"

    # Weight specific
    if "体重" in t:
        if "增加" in direction or "高" in direction:
            gain_kw = ("体重增加", "水肿", "食欲增加", "超重", "肥胖")
            if any(k in combined_text for k in gain_kw):
                return True, "体重增加相关AE/MH可解释"
        if "降低" in direction or "低" in direction:
            loss_kw = ("体重降低", "体重下降", "消瘦", "食欲减退", "厌食", "营养不良", "腹泻", "呕吐", "胃炎", "消化不良")
            if any(k in combined_text for k in loss_kw):
                return True, "体重降低或胃肠/营养相关AE/MH可解释"

    # General relevance rules
    rules = [
        (("白细胞", "WBC"), ("白细胞", "中性粒", "感染", "咽喉炎", "上呼吸道感染", "炎症", "咽炎", "肺炎", "鼻窦炎"), "感染/炎症可解释白细胞波动"),
        (("中性粒细胞", "NEUT"), ("中性粒", "感染", "炎症", "白细胞"), "感染/炎症可解释中性粒细胞波动"),
        (("淋巴细胞", "LYMPH"), ("淋巴细胞", "感染", "炎症", "病毒感染"), "感染/炎症可解释淋巴细胞波动"),
        (("嗜酸性粒细胞", "EO"), ("嗜酸", "过敏", "皮炎", "鼻息肉", "鼻窦炎", "本研究疾病"), "过敏/炎症可解释嗜酸性粒细胞波动"),
        (("嗜碱性粒细胞", "BASO"), ("嗜碱", "过敏", "皮炎", "炎症"), "过敏/炎症可解释嗜碱性粒细胞波动"),
        (("血小板", "PLT"), ("血小板", "PLT"), "血小板相关AE/MH可解释"),
        (("血红蛋白", "Hb", "血细胞比容", "HCT", "红细胞", "RBC"), ("贫血", "血红蛋白", "红细胞"), "贫血相关AE/MH可解释"),
        (("总胆固醇", "CHOL", "甘油三酯", "TG"), ("高脂血症", "血脂", "胆固醇", "甘油三酯", "脂肪肝"), "血脂异常相关AE/MH可解释"),
        (("ALT", "AST", "转氨酶"), ("转氨酶", "肝功能", "脂肪肝", "肝炎", "肝"), "肝功能相关AE/MH可解释"),
        (("碱性磷酸酶", "ALP"), ("碱性磷酸酶", "ALP", "肝", "胆"), "ALP/肝胆相关AE/MH可解释"),
        (("总胆红素", "TBiL", "直接胆红素", "DBiL", "间接胆红素", "IBil"), ("胆红素", "黄疸", "肝", "胆"), "胆红素/肝胆相关AE/MH可解释"),
        (("白蛋白", "ALB"), ("白蛋白", "低蛋白", "营养", "肝"), "白蛋白/营养相关AE/MH可解释"),
        (("肌酐", "Cr"), ("肌酐", "肾", "肾功能", "肾损伤"), "肌酐/肾功能相关AE/MH可解释"),
        (("尿素氮", "BUN", "血尿素", "Urea"), ("尿素", "肾", "肾功能"), "尿素/肾功能相关AE/MH可解释"),
        (("血葡萄糖", "GLU"), ("血糖", "葡萄糖", "糖尿病", "高血糖", "低血糖"), "血糖相关AE/MH可解释"),
        (("血清钾", "K+"), ("钾", "低钾", "高钾", "电解质"), "电解质相关AE/MH可解释"),
        (("血清钠", "Na+"), ("钠", "低钠", "高钠", "电解质"), "电解质相关AE/MH可解释"),
        (("血清钙", "Ca"), ("钙", "低钙", "高钙", "电解质"), "电解质相关AE/MH可解释"),
        (("血清镁", "Mg"), ("镁", "低镁", "高镁", "电解质"), "电解质相关AE/MH可解释"),
        (("总蛋白", "TP"), ("蛋白", "营养", "肝"), "蛋白/营养相关AE/MH可解释"),
        (("收缩压",), ("高血压", "低血压", "血压"), "血压相关AE/MH可解释"),
        (("舒张压",), ("高血压", "低血压", "血压"), "血压相关AE/MH可解释"),
        (("QTcF", "QT间期", "QT"), ("QT", "心电", "心律", "心动过速", "心动过缓", "心律失常"), "心电/心律相关AE/MH可解释"),
        (("心率",), ("心率", "心电", "心律", "心动过速", "心动过缓"), "心率相关AE/MH可解释"),
    ]

    for needles, keys, reason in rules:
        if any(n.lower() in t.lower() for n in needles):
            if any(k.lower() in combined_text.lower() for k in keys):
                return True, reason

    return False, ""


def parse_ae_mh_from_desc(desc):
    """Extract AE/MH numbers from description text."""
    if not desc:
        return []
    refs = []
    for m in re.finditer(r"(AE|MH)\s*(\d+)", desc, re.IGNORECASE):
        refs.append((m.group(1).upper(), int(m.group(2))))
    return refs


def assess_risk(parsed_data, ctcae_index, non_ctcae_threshold=20):
    """
    Main risk assessment function.
    Returns categorized results.
    """
    lab_records = parsed_data["lab_records"]
    ae_records = parsed_data.get("ae_records", [])
    mh_records = parsed_data.get("mh_records", [])

    # Build AE/MH lookup
    ae_lookup = defaultdict(list)
    for r in ae_records:
        ae_lookup[r["subj"]].append(r)
    mh_lookup = defaultdict(list)
    for r in mh_records:
        mh_lookup[r["subj"]].append(r)

    # Group records by (center, subject, test)
    by_subject_test = defaultdict(list)
    for rec in lab_records:
        key = (rec.get("center", ""), rec["subj"], rec["test"])
        by_subject_test[key].append(rec)

    # Sort each group by visit
    for key in by_subject_test:
        by_subject_test[key].sort(
            key=lambda r: (get_visit_rank(r.get("visit", "")), r.get("date") or "")
        )

    results = defaultdict(lambda: defaultdict(list))

    for (center, subj, test), records in by_subject_test.items():
        # Find baseline
        baselines = [r for r in records if is_baseline_visit(r.get("visit", ""))]
        if not baselines:
            baselines = [r for r in records if is_screening_visit(r.get("visit", ""))]
        if not baselines:
            baselines = [records[0]]

        baseline = baselines[-1]  # Use last baseline row
        base_result = baseline.get("result")
        base_visit = baseline.get("visit", "")

        # Check if this test has CTCAE
        has_ctcae = ctcae_index.get(test, {}).get("has_grades", False)

        # Identify abnormal post-baseline records
        abnormal_posts = []
        all_visits = []
        for rec in records:
            lo = rec.get("lo")
            hi = rec.get("hi")
            result = rec.get("result")
            phase = classify_visit(rec.get("visit", ""))

            is_ab, direction = is_abnormal(result, lo, hi)
            grade, grade_desc = get_ctcae_grade(test, result, lo, hi, base_result, ctcae_index)

            # For non-CTCAE tests, check percentage change
            if not has_ctcae and phase == "post" and base_result:
                pct = None
                if base_result != 0 and isinstance(result, (int, float)):
                    pct = (result - base_result) / base_result * 100
                # Check if significant enough
                if pct and abs(pct) >= non_ctcae_threshold:
                    is_ab = True
                    direction = "高" if pct > 0 else "低"
                elif not is_ab:
                    is_ab = False
                    direction = None

            # CTCAE upgrade check
            if baseline in baselines:
                base_grade, _ = get_ctcae_grade(test, base_result, baseline.get("lo"),
                                                 baseline.get("hi"), base_result, ctcae_index)
            else:
                base_grade = None

            # Weight matching
            if test in ("体重",) and grade and grade >= 1:
                is_ab = True

            # Determine if CTCAE upgraded from baseline
            ctcae_upgraded = bool(grade and base_grade and grade > base_grade)

            # Rule 1: Baseline abnormal → baseline normal → NOT a risk item
            bl_is_ab, _ = is_abnormal(base_result, baseline.get("lo"), baseline.get("hi"))
            if bl_is_ab and not is_ab and not ctcae_upgraded:
                # Baseline was abnormal, post is normal, no CTCAE upgrade → skip
                is_ab = False

            # Combine desc from multiple sources: desc, desc_link, note
            combined_desc = rec.get("desc", "") or rec.get("desc_link", "") or rec.get("note", "")

            rec_info = {
                "visit": rec.get("visit", ""),
                "visit_point": rec.get("visit_point", ""),
                "date": rec.get("date"),
                "result": result,
                "lo": lo,
                "hi": hi,
                "unit": rec.get("unit", ""),
                "is_ab": is_ab,
                "direction": direction,
                "cs": rec.get("cs", ""),
                "desc": combined_desc,
                "ae_ref": rec.get("ae_ref", ""),
                "mh_ref": rec.get("mh_ref", ""),
                "grade": grade,
                "grade_desc": grade_desc,
                "ctcae_upgraded": ctcae_upgraded,
                "phase": phase,
            }

            # Calculate percentage change from baseline
            if phase == "post" and base_result and isinstance(result, (int, float)) and base_result != 0:
                rec_info["pct_change"] = (result - base_result) / base_result * 100
            else:
                rec_info["pct_change"] = None

            all_visits.append(rec_info)
            if phase == "post" and (is_ab or ctcae_upgraded):
                abnormal_posts.append(rec_info)

        if not abnormal_posts:
            continue

        # Collect AE/MH from description references
        linked_ae = []
        linked_mh = []
        for rec_info in abnormal_posts:
            desc = rec_info.get("desc", "")
            for kind, num in parse_ae_mh_from_desc(desc):
                if kind == "AE":
                    for ae in ae_lookup.get(subj, []):
                        linked_ae.append(ae)
                elif kind == "MH":
                    for mh in mh_lookup.get(subj, []):
                        linked_mh.append(mh)

        # Also check AE/MH reference columns
        for rec_info in abnormal_posts:
            ae_ref = rec_info.get("ae_ref", "")
            mh_ref = rec_info.get("mh_ref", "")
            if ae_ref:
                for ae in ae_lookup.get(subj, []):
                    linked_ae.append(ae)
            if mh_ref:
                for mh in mh_lookup.get(subj, []):
                    linked_mh.append(mh)

        # Deduplicate
        def dedup(items, key="term"):
            seen = set()
            result = []
            for item in items:
                k = item.get(key, "")
                if k and k not in seen:
                    seen.add(k)
                    result.append(item)
            return result

        linked_ae = dedup(linked_ae)
        linked_mh = dedup(linked_mh)

        # Determine AE/MH relevance
        ae_terms = [a.get("term", "") for a in linked_ae]
        mh_terms = [m.get("term", "") for m in linked_mh]
        all_terms = ae_terms + mh_terms

        # Determine direction
        direction = ""
        for p in abnormal_posts:
            if p.get("direction"):
                direction = p["direction"]
                break

        relevant, reason = term_relevance(test, all_terms, direction)

        # Count risk factors - proper CS classification
        cs_normal_count = sum(1 for p in abnormal_posts
                              if norm(p.get("cs", "")) in ("正常", "Normal", "normal"))
        ncs_count = sum(1 for p in abnormal_posts
                        if "无临床意义" in p.get("cs", "") or "NCS" in norm(p.get("cs", "")))
        cs_count = sum(1 for p in abnormal_posts
                       if "有临床意义" in p.get("cs", "") or "CS" in norm(p.get("cs", ""))
                       if "无临床意义" not in p.get("cs", "") and "NCS" not in norm(p.get("cs", "")))
        no_eval_count = sum(1 for p in abnormal_posts
                           if not p.get("cs") or norm(p.get("cs")) in ("", "-", "nan"))

        grade_count = sum(1 for p in abnormal_posts if p.get("grade"))
        ctcae_upgrade_count = sum(1 for p in abnormal_posts if p.get("ctcae_upgraded"))
        pct10_count = sum(1 for p in abnormal_posts
                         if p.get("pct_change") and abs(p["pct_change"]) >= 10)

        # Determine risk level using refined logic (v2 - user's spec)
        # Risk classification per spec:
        # (A) CS="异常无临床意义"(NCS) or unevaluated → HIGH risk
        # (B) Mixed NCS + CS → MEDIUM risk
        # (C) All CS with inconsistent explanation → MEDIUM risk
        # (D) CS with consistent AE/MH explanation → LOW risk
        # (E) CS="正常"(technically abnormal but investigator says normal) → LOW risk

        total_ab = len(abnormal_posts)

        # If all abnormal posts have CS="正常": technically outside range but investigator normal
        if cs_normal_count == total_ab:
            risk = "low"
        # Mixed: some CS normal, some NCS
        elif cs_normal_count > 0 and ncs_count > 0:
            risk = "medium"
        # NCS or unevaluated → HIGH (per spec)
        elif ncs_count > 0 or no_eval_count > 0:
            if relevant and (linked_ae or linked_mh):
                risk = "medium"  # AE/MH linked but still NCS
            else:
                risk = "high"
        # All CS with consistent AE/MH explanation
        elif cs_count == total_ab and relevant and linked_ae:
            risk = "low"
        # All CS but no AE/MH or not relevant
        elif cs_count == total_ab:
            risk = "medium"
        # Mixed CS + NCS
        elif cs_count > 0 and ncs_count > 0:
            risk = "medium"
        # Default
        else:
            risk = "medium"

        # Build AE/MH info for card
        ae_info = []
        for ae in linked_ae:
            ae_info.append({
                "term": ae.get("term", ""),
                "severity": ae.get("severity", ""),
                "relationship": ae.get("relationship", ""),
                "start": ae.get("start"),
                "outcome": ae.get("outcome"),
            })
        mh_info = []
        for mh in linked_mh:
            mh_info.append({
                "term": mh.get("term", ""),
                "start": mh.get("start"),
                "ongoing": mh.get("ongoing", ""),
            })

        # Generate explain evaluation and response
        if not linked_ae and not linked_mh:
            explain_eval = "未见可核对的AE/MH解释"
        elif relevant:
            explain_eval = reason
        else:
            explain_eval = "已关联AE/MH但与本指标异常关联不足"

        # Risk items
        item = {
            "site": center,
            "subj": subj,
            "test": test,
            "has_ctcae": has_ctcae,
            "ctcae_term": ctcae_index.get(test, {}).get("ctcae_term", ""),
            "risk": risk,
            "bl": {
                "visit": base_visit,
                "date": baseline.get("date"),
                "result": base_result,
                "unit": baseline.get("unit", ""),
                "lo": baseline.get("lo"),
                "hi": baseline.get("hi"),
            },
            "bl_status": "基线正常" if not is_abnormal(base_result,
                baseline.get("lo"), baseline.get("hi"))[0] else "基线异常",
            "bl_dir": is_abnormal(base_result, baseline.get("lo"), baseline.get("hi"))[1] or "",
            "post": abnormal_posts,
            "all_visits": all_visits,
            "ae": ae_info,
            "mh": mh_info,
            "cm": [],
            "excluded_ae": [],
            "excluded_mh": [],
            "new_ab_count": len(abnormal_posts),
            "grade_count": grade_count,
            "clin_sig_count": cs_count,
            "cs_normal_count": cs_normal_count,
            "ncs_count": ncs_count,
            "no_eval_count": no_eval_count,
            "ctcae_upgrade_count": ctcae_upgrade_count,
            "baseline_change_10pct_count": pct10_count,
            "explain_eval": explain_eval,
        }
        results[center][test].append(item)

    # Sort by risk within each test
    for center in results:
        for test in results[center]:
            results[center][test].sort(
                key=lambda x: {"high": 0, "medium": 1, "low": 2}.get(x["risk"], 9)
            )

    return results


def run_assessment(parsed_data_path, ctcae_index_path=None, output_path=None,
                   non_ctcae_threshold=20):
    """Run full assessment pipeline."""
    parsed_data = json.loads(Path(parsed_data_path).read_text(encoding="utf-8"))

    if ctcae_index_path:
        ctcae_index = json.loads(Path(ctcae_index_path).read_text(encoding="utf-8"))
        ctcae_mapping = ctcae_index.get("mapping", {})
    else:
        from parse_ctcae import CTCAE_V5_MAPPING
        ctcae_mapping = CTCAE_V5_MAPPING

    from parse_ctcae import build_ctcae_index
    ctcae_index = build_ctcae_index(parsed_data["test_names"], ctcae_mapping)

    results = assess_risk(parsed_data, ctcae_index, non_ctcae_threshold)

    # Summary stats
    stats = {}
    for center, tests in results.items():
        all_items = [it for items in tests.values() for it in items]
        stats[center] = {
            "high": sum(1 for it in all_items if it["risk"] == "high"),
            "medium": sum(1 for it in all_items if it["risk"] == "medium"),
            "low": sum(1 for it in all_items if it["risk"] == "low"),
            "subjects": len(set(it["subj"] for it in all_items)),
            "total": len(all_items),
        }

    output = {
        "results": results,
        "stats": stats,
        "ctcae_index": ctcae_index,
        "non_ctcae_threshold": non_ctcae_threshold,
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "test_names": parsed_data["test_names"],
        "test_sources": parsed_data["test_sources"],
    }

    if output_path:
        Path(output_path).write_text(
            json.dumps(output, ensure_ascii=False, indent=2),
            encoding="utf-8"
        )
        print(f"Assessment saved to {output_path}")

    # Print summary
    for center, s in stats.items():
        print(f"  {center}: 高{s['high']} / 中{s['medium']} / 低{s['low']} / 受试者{s['subjects']} / 总{s['total']}")

    return output


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python assess_risk.py <parsed_data.json> [ctcae_index.json] [output.json] [threshold%]")
        sys.exit(1)

    threshold = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    run_assessment(
        sys.argv[1],
        sys.argv[2] if len(sys.argv) > 2 else None,
        sys.argv[3] if len(sys.argv) > 3 else None,
        threshold,
    )
