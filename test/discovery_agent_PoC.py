import os

# 必须在导入 langchain 之前设置，完全禁用 LangSmith 追踪
os.environ["LANGCHAIN_TRACING_V2"] = "false"
os.environ["LANGCHAIN_TRACING"] = "false"
os.environ["LANGCHAIN_ENDPOINT"] = ""
os.environ["LANGCHAIN_API_KEY"] = ""
os.environ["LANGCHAIN_PROJECT"] = ""

import json
import re
import time
from typing import Dict, List
from xml.etree import ElementTree as ET

import pandas as pd
from Bio import Entrez
from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.vectorstores.base import VectorStoreRetriever
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_openai import ChatOpenAI
from pydantic import BaseModel, Field

# api_key = "sk-KaZVAPnsPr2oVbLq17511e02E979454bBd43E0B07b18344f"
# base_url = "https://api.pumpkinaigc.online/v1"
# llm = ChatOpenAI(api_key=api_key, base_url=base_url, model="gpt-4.1")


# 1.抽取PICO元素（使用结构化输出）
# class PICOModel(BaseModel):
#     """PICO模型：临床研究问题框架"""

#     P: str = Field(default="", description="Population - 研究人群")
#     I: str = Field(default="", description="Intervention - 干预措施")
#     C: str = Field(default="", description="Comparison - 对照组")
#     O: str = Field(default="", description="Outcome - 结局指标")


# def extract_pico(text):
#     structured_llm = llm.with_structured_output(PICOModel)

#     prompt = ChatPromptTemplate.from_messages(
#         [
#             (
#                 "system",
#                 "你是临床研究专家。从文本中提取PICO要素（Population、Intervention、Comparison、Outcome）。如果某个元素不存在，留空。",
#             ),
#             ("user", "文本知识：{text}"),
#         ]
#     )

#     chain = prompt | structured_llm
#     result = chain.invoke({"text": text})

#     return result.model_dump()


# 2.基于Langchain的PubMedRetriever检索临床研究领域文献,抽取摘要中的PICO，将其写入PICO知识库
# def search(keyword, max_res_count):
#     while True:
#         try:
#             time.sleep(0.2)
#             handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_res_count)
#             record = Entrez.read(handle)
#             return record["IdList"]
#         except Exception as e:
#             print("search error:", e)
#             time.sleep(2)


# def fetch(pmid):
#     while True:
#         try:
#             time.sleep(0.2)
#             handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml")
#             xml_data = handle.read()

#             root = ET.fromstring(xml_data)
#             abstract_text = ""
#             for abstract in root.findall(".//Abstract"):
#                 for abstract_text_elem in abstract.findall("AbstractText"):
#                     abstract_text += abstract_text_elem.text
#             clean_abstract = re.sub(r"[\n\r\t]", "", abstract_text)
#             return clean_abstract.strip() or "No abstract available."

#         except Exception as e:
#             print("fetch error:", e)
#             time.sleep(2)


# pmids = search("clinical trials", 10)
# abstracts = []
# for pmid in pmids:
#     abs_text = fetch(pmid)
#     abstracts.append({pmid: abs_text})
# print("成功获取摘要数量:", len(abstracts))
# print(abstracts[:1])

# knowledge_base = []
# for item in abstracts:
#     for pmid, abstract in item.items():
#         pico_entry = {
#             **extract_pico(abstract),
#             "id": len(knowledge_base) + 1,
#             "pmid": pmid,
#             "abstract": abstract,
#         }
#         knowledge_base.append(pico_entry)

# df = pd.DataFrame(knowledge_base)
# df.head()
# df.to_csv("pico_knowledge_base.csv", index=False)

# with open("pico_knowledge_base.json", "w", encoding="utf-8") as f:
#     for item in knowledge_base:
#         f.write(json.dumps(item, ensure_ascii=False) + "\n")

# 3.向量化PICO知识库
# model = SentenceTransformer('BAAI/bge-small-en-v1.5')
# embedding_model = HuggingFaceEmbeddings(
#     model_name="BAAI/bge-small-en-v1.5",
#     model_kwargs={"device": "cpu"},
#     encode_kwargs={"normalize_embeddings": True},
# )
# docs = []
# for record in knowledge_base:
#     meta = {"pmid": record.get("pmid")}
#     text = f"""
#                 PMID: {record.get("pmid")};
#                 P: {record.get("P")};
#                 I: {record.get("I")};
#                 C: {record.get("C")};
#                 O: {record.get("O")};
#                 """
#     clean_text = re.sub(r"[\n\r\t]", "", text).strip()
#     print(clean_text)
#     docs.append(Document(page_content=clean_text, metadata=meta))
faiss_index = FAISS.from_documents(docs, embedding_model)
retriever = faiss_index.as_retriever(search_kwargs={"k": 2})


# 4.抽取用户问题的PICO（以不完整PICO为例）,转换为自然语言文本，在PICO知识库中进行相似性检索，找出相关的PICOs
def build_query_text(partial_pico: Dict[str, str]) -> str:
    parts = []
    for slot in ["P", "I", "C", "O"]:
        val = partial_pico.get(slot) or partial_pico.get(slot.lower()) or ""
        if val:
            parts.append(f"{slot}: {val}")
    return "\n".join(parts)


def get_relevant_docs(query_text: str, retriever: VectorStoreRetriever) -> List[Document]:
    qtxt = build_query_text(query_text)
    return retriever.invoke(qtxt)


# 5.利用LLM推理填充用户问题的PICO
complement_prompt = """
你是临床研究领域专家。

现在有一个用户给出的不完整的问题模型JSON（PICO模型：population、intervention，comparison，outcome）：
{user_query_pico}

在知识库中检索到与用户问题相关的一个或多个证据问题模型列表：
{evidence_pico_docs}

请完成以下任务：

1. 列出“用户问题模型文本”与“证据问题模型文本”在 P、I、C、O 四个维度上的相似性与差异性。
2. 推断“用户问题模型文本”缺失或可能延伸的研究空白点。
3. 基于推断的空白点与“证据问题模型文本”，填充“用户问题模型文本”，生成5个完整的PICO模型，要求兼具创新性、可行性与科学性。

输出JSON格式，用中文回答。
"""

user_query = "I want to study weight outcomes in Type 2 diabetes. Any gaps?"
user_query_pico = extract_pico(user_query)
evidence_pico_docs = get_relevant_docs(user_query_pico, retriever)

# 使用 prompt 和 LLM 生成答案
formatted_prompt = complement_prompt.format(
    user_query_pico=user_query_pico, evidence_pico_docs=evidence_pico_docs
)
result = llm.invoke([{"role": "user", "content": formatted_prompt}])

print("Answer:\n", result.content)
print("Source docs:\n", [d.metadata for d in docs])


# from langchain_community.retrievers import PubMedRetriever
# import numpy as np

# import os
# import urllib.request
# os.environ["HTTPS_PROXY"] = "http://127.0.0.1:7890"
# try:
#     response = urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi", timeout=10)
#     print("✅ PubMed 接口可访问！")
# except Exception as e:
#     print(f"❌ 接口访问失败：{str(e)}")

# def retrievePubMedByLangchain(keyword, max_res_count):
#     retriever = PubMedRetriever(top_k=max_res_count)
#     try:
#         docs = retriever.invoke(keyword)
#         print(f"成功获取 {len(docs)} 篇文献")
#     except Exception as e:
#         print(f"失败获取：{str(e)}")
#     abstracts = []
#     for doc in docs:
#         abstract = doc.page_content
#         pmid = doc.metadata.get("uid", "")
#         abstracts.append({pmid: abstract})
#     return abstracts

# abstracts = retrievePubMed("clinical trials", 1)
# list(abstracts.items())[:2]
