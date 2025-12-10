"""
测试 extract_research_type 函数
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from question_model_extractor import QuestionModelExtractor

# 创建提取器实例
extractor = QuestionModelExtractor()

# 测试案例1: 随机对照试验
abstract1 = {
    "title": "Effect of Metformin on HbA1c in Type 2 Diabetes: A Randomized Controlled Trial",
    "abstract": """
    Background: This study aimed to evaluate the efficacy of metformin in controlling blood glucose levels.
    Methods: A total of 200 patients with type 2 diabetes were randomly assigned to receive either 
    metformin 1000mg twice daily or placebo for 12 weeks. The primary outcome was change in HbA1c.
    Results: The metformin group showed a mean HbA1c reduction of 1.5% compared to 0.3% in the placebo group (p<0.001).
    Conclusion: Metformin significantly improves glycemic control in type 2 diabetes patients.
    """
}

# 测试案例2: 队列研究
abstract2 = {
    "title": "Long-term cardiovascular outcomes in diabetic patients: A 10-year cohort study",
    "abstract": """
    Objective: To investigate the incidence of cardiovascular events in patients with diabetes.
    Methods: We followed 1500 diabetic patients from 2010 to 2020, recording all cardiovascular events.
    The cohort was stratified by baseline HbA1c levels. Cox regression was used to analyze risk factors.
    Results: Higher baseline HbA1c was associated with increased cardiovascular risk (HR 1.8, 95% CI 1.4-2.3).
    """
}

# 测试案例3: 病例对照研究
abstract3 = {
    "title": "Risk factors for diabetic retinopathy: A case-control study",
    "abstract": """
    Purpose: To identify risk factors associated with diabetic retinopathy.
    Methods: We compared 150 diabetic patients with retinopathy (cases) to 150 diabetic patients 
    without retinopathy (controls). Logistic regression was performed to identify risk factors.
    Results: Duration of diabetes, poor glycemic control, and hypertension were significant risk factors.
    """
}

# 测试案例4: 交叉试验
abstract4 = {
    "title": "Comparison of two insulin regimens in type 1 diabetes: A crossover trial",
    "abstract": """
    Aim: To compare the efficacy of two different insulin regimens.
    Design: In this crossover trial, 50 type 1 diabetes patients received regimen A for 3 months,
    then switched to regimen B for another 3 months after a 2-week washout period.
    Results: Both regimens showed similar efficacy in glucose control.
    """
}

print("=" * 80)
print("测试 extract_research_type 函数")
print("=" * 80)

test_cases = [
    ("随机对照试验", abstract1),
    ("队列研究", abstract2),
    ("病例对照研究", abstract3),
    ("交叉试验", abstract4)
]

for name, abstract in test_cases:
    print(f"\n{'='*80}")
    print(f"测试案例: {name}")
    print(f"{'='*80}")
    print(f"标题: {abstract['title']}")
    print(f"\n提取的研究类型:")
    
    research_type = extractor.extract_research_type(abstract)
    print(f"✓ {research_type}")
    print()

print("=" * 80)
print("完整流程示例")
print("=" * 80)

# 完整流程：任务类型 → 问题模型 → 研究设计类型
abstract = abstract1

print("\n1. 判断任务类型:")
task_type = extractor.determine_task_type(abstract)
print(f"   任务类型: {task_type}")

print("\n2. 提取研究设计类型:")
research_type = extractor.extract_research_type(abstract)
print(f"   研究设计类型: {research_type}")

print("\n3. 提取问题模型 (需要数据库连接):")
print(f"   根据任务类型 '{task_type}' 应该使用的问题模型:")
from constants import TASK_TO_MODEL_MAP
model_type = TASK_TO_MODEL_MAP.get(task_type, "Unknown")
print(f"   问题模型: {model_type}")

print("\n" + "=" * 80)
print("测试完成!")
print("=" * 80)
