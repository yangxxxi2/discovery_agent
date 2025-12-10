DETERMINE_TASK_TYPE_PROMPT = """
你是一名循证医学专家，请分析以下医学文献摘要与题目，判断其任务类型。

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

文献信息：
标题: {title}
摘要: {abstract}

请从以上选项中选择最符合的研究设计类型，只返回类型名称（英文），不要包含其他解释。
如果无法判断，返回 "Unknown"。
"""

EXTRACT_QUESTION_MODEL_PROMPT = """
请先理解{model_type}模型的相关信息，从以下由医学文献摘要与题目组成的文本中提取该模型的各个组件。

{model_type}模型的相关信息：
1. 该模型主要适用于临床研究领域中的{model_description}
2. 适用于该模型的文本例子：{model_examples}
3. 该模型包含以下组件：
{components_info}

医学文献摘要与题目组成的文本：
{text}

请以JSON格式返回提取结果，如果某个组件在文本中没有明确提及，请填写"None"。

示例格式：
{example_format}

请直接返回JSON，不要包含其他文字。
"""

EXTRACT_RESEARCH_TYPE_PROMPT = """
你是一名循证医学专家，请分析以下医学文献摘要与题目，判断其研究设计类型。

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

文献信息：
标题: {title}
摘要: {abstract}

请从以上选项中选择最符合的研究设计类型，只返回类型名称（英文），不要包含其他解释。
如果无法判断，返回 "Unknown"。
"""