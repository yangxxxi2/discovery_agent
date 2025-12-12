DETERMINE_TASK_TYPE_PROMPT = """
你是一名循证医学专家，请分析以下文本，判断其任务类型。
这个文本有可能是医学文献的标题和摘要，也有可能是人类对研究问题的描述。

任务类型选项：
1. Interventional Task with Time Dimension - 干预性研究，包含时间维度
2. Interventional Task without Time Dimension - 干预性研究，不包含时间维度
3. Observational Task with Comparison - 观察性研究，包含对照组
4. Observational Task without Comparison - 观察性研究，不包含对照组
5. Prevalence / Incidence Task - 流行率/发病率研究
6. Diagnostic Test Accuracy Task - 诊断测试准确性研究
7. Qualitative Task - 质性研究
8. Health Practice Task - 卫生实践研究
9. Health Service Task - 卫生服务研究

文本：
{text}

请从以上选项中选择最符合的研究设计类型，只返回类型名称（英文），不要包含其他解释。
如果无法判断，返回 "Unknown"。
"""

EXTRACT_QUESTION_MODEL_PROMPT = """
请先理解{model_type}模型的相关信息，从以下文本中提取该模型的各个组件。
这个文本有可能是医学文献的标题和摘要，也有可能是人类对研究问题的描述。

{model_type}模型的相关信息：
1. 该模型主要适用于临床研究领域中的{model_description}
2. 适用于该模型的文本例子：{model_examples}
3. 该模型包含以下组件：
{components_info}

文本：
{text}

请以JSON格式返回提取结果，如果某个组件在文本中没有明确提及，请填写"Unknown"。

示例格式：
{example_format}

请直接返回JSON，不要包含其他文字。
"""

EXTRACT_RESEARCH_TYPE_PROMPT = """
你是一名循证医学专家，请分析以下文本，判断其研究设计类型。
这个文本有可能是医学文献的标题和摘要，也有可能是人类对研究问题的描述。

研究设计类型选项：
1. Randomized Controlled Trial - 随机对照试验
2. Non Randomized Controlled Trial - 非随机对照试验
3. Single Arm Trial - 单臂试验
4. Three Arm Trial - 三臂试验
5. Crossover Trial - 交叉试验
6. Factorial Trial - 析因试验
7. Cohort Study - 队列研究
8. Case-Control Study - 病例对照研究
9. Cross-Sectional Study - 横断面研究
10. Before-and-After Study - 前后对照研究
11. Ecological Study - 生态学研究
12. Focus Group Discussion - 焦点小组讨论
13. In-depth Interview - 深度访谈

文本：
{text}

请从以上选项中选择最符合的研究设计类型，只返回类型名称（英文），不要包含其他解释。
如果无法判断，返回 "Unknown"。
"""

EXTEND_FILLED_SLOT_PROMPT = """
你是一名循证医专家。请为以下研究问题组件生成语义相似的变体表述。

研究框架：{pending_model.get('model_type', 'Unknown')}
研究类型：{pending_model.get('research_type', 'Unknown')}

已填充的组件：
{json.dumps(filled_slots, ensure_ascii=False, indent=2)}

任务要求：
1. 为每个组件生成2-3个语义相似但表述不同的变体
2. 变体应该：
   - 保持核心含义不变
   - 使用不同的医学术语或表达方式
   - 涵盖相关的概念范围
3. 有助于扩大文献检索的覆盖面

返回JSON格式：
{{
  "组件名": ["变体1", "变体2", "变体3"],
  ...
}}

示例：
输入: {{"P": "青少年抑郁症"}}
输出: {{"P": ["青春期抑郁障碍", "未成年人抑郁症状", "少年期情绪障碍"]}}
"""

ANALOGICAL_GENERATION_PROMPT = """
你是一名循证医学专家。基于以下信息，通过举一反三生成创新性的研究问题。

用户的原始问题（部分信息）：
{json.dumps(pending_model, ensure_ascii=False, indent=2)}

已有的相关研究（供参考，但不要重复）：
{json.dumps([{{
    'title': e['title'],
    'research_type': e['research_type'],
    'components': e['components'],
    'similarity': f"{e['similarity_score']:.2f}"
}} for e in evidence[:5]], ensure_ascii=False, indent=2)}

任务要求：
1. 分析已有研究的模式和特点
2. 识别研究gap（哪些方向还未被充分探索）
3. 基于用户问题，生成3-5个创新性的完整研究问题
4. 这些问题应该：
   - 与已有研究相关但不重复
   - 有临床意义和研究价值
   - 可行且符合伦理
   - 填补现有研究的空白

请以JSON数组格式返回，每个问题包含所有必要的组件。

示例格式：
[
  {{
    "model_type": "PICO",
    "P": "具体人群",
    "I": "具体干预",
    "C": "具体对照",
    "O": "具体结局",
    "innovation_rationale": "创新点说明"
  }}
]
"""

EVALUATE_CANDIDATE_MODELS_PROMPT = """
你是一名资深的循证医学评审专家。请评估以下研究问题的质量。

候选研究问题：
{json.dumps(candidate_models, ensure_ascii=False, indent=2)}

评估标准（1-10分）：
1. 创新性 (Novelty): 是否提出了新的研究视角？
2. 可行性 (Feasibility): 研究是否可以实际执行？
3. 临床意义 (Clinical Significance): 对临床实践的潜在价值？
4. 完整性 (Completeness): 问题表述是否清晰完整？
5. 伦理性 (Ethical): 是否符合研究伦理？

请为每个问题评分并排序，返回JSON格式：
[
  {{
    "model": {{...原始模型...}},
    "scores": {{
      "novelty": 分数,
      "feasibility": 分数,
      "clinical_significance": 分数,
      "completeness": 分数,
      "ethical": 分数
    }},
    "total_score": 总分,
    "rationale": "评分理由"
  }}
]

按total_score从高到低排序。


"""