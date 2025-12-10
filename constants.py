#任务类型
TASK_TYPE_INTEVENTION_TIME = 'Interventional Task with Time Dimension'
TASK_TYPE_INTEVENTION = 'Interventional Task without Time Dimension'
TASK_TYPE_OBSERVATION_COMPARISON = 'Observational Task with Comparison'
TASK_TYPE_OBSERVATION = 'Observational Task without Comparison'
TASK_TYPE_PREVALENCE = 'Prevalence / Incidence Task'
TASK_TYPE_DIAGNOSIS = 'Diagnostic Test Accuracy Task'
TASK_TYPE_QUALITATION = 'Qualitative Task'
TASK_TYPE_HEALTH_PRACTICE = 'Health Practice Task'
TASK_TYPE_HEALTH_SERVICE = 'Health Service Task'

#问题模型
QUESTION_MODEL_PICO = 'PICO'
QUESTION_MODEL_PICOT = 'PICOT'
QUESTION_MODEL_PECO = 'PECO'
QUESTION_MODEL_PEO = 'PEO'
QUESTION_MODEL_CoCoPop = 'CoCoPop'
QUESTION_MODEL_PIRD = 'PIRD'
QUESTION_MODEL_SPIDER = 'SPIDER'
QUESTION_MODEL_SPICE = 'SPICE'
QUESTION_MODEL_ECLIPSe = 'ECLIPSe'

#任务类型到问题模型映射
TASK_TO_MODEL_MAP = {
            TASK_TYPE_INTEVENTION_TIME: QUESTION_MODEL_PICOT,
            TASK_TYPE_INTEVENTION: QUESTION_MODEL_PICO,
            TASK_TYPE_OBSERVATION_COMPARISON: QUESTION_MODEL_PECO,
            TASK_TYPE_OBSERVATION: QUESTION_MODEL_PEO,
            TASK_TYPE_PREVALENCE: QUESTION_MODEL_CoCoPop,
            TASK_TYPE_DIAGNOSIS: QUESTION_MODEL_PIRD,
            TASK_TYPE_QUALITATION: QUESTION_MODEL_SPIDER,
            TASK_TYPE_HEALTH_PRACTICE: QUESTION_MODEL_SPICE,
            TASK_TYPE_HEALTH_SERVICE: QUESTION_MODEL_ECLIPSe,
        }

#研究设计类型
RESEARCH_TYPE_COHORT = 'Cohort Study'
RESEARCH_TYPE_CASE_CONTROL = 'Case-Control Study'
RESEARCH_TYPE_CROSS_SECTIONAL = 'Cross-Sectional Study'
RESEARCH_TYPE_RCT = 'Randomized Controlled Trial'
RESEARCH_TYPE_NRCT = 'Non Randomized Controlled Trial'
RESEARCH_TYPE_SAT = 'Single Arm Trial'
RESEARCH_TYPE_TAT = 'Three Arm Trial'
RESEARCH_TYPE_CROSS = 'Crossover Trial'
RESEARCH_TYPE_FACTORIAL = 'Factorial Trial'
RESEARCH_TYPE_FGD = 'Focus Group Discussion'
RESEARCH_TYPE_IN_DEPTH = 'In-depth Interview'
RESEARCH_TYPE_BEFORE_AFTER = 'Before-and-After Study'
RESEARCH_TYPE_ECOLOGICAL = 'Ecological Study'

#检索关键字
KEYWORD_CLINICAL_RESEARCH = "clinical research OR clinical study OR clinical trial OR randomized controlled trial OR cohort study OR case-control study OR observational study OR diagnostic study"