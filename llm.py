from langchain_openai import ChatOpenAI

api_key = "sk-KaZVAPnsPr2oVbLq17511e02E979454bBd43E0B07b18344f"
base_url = "https://api.pumpkinaigc.online/v1"
llm = ChatOpenAI(api_key=api_key, base_url=base_url, model="gpt-4.1", temperature=0)
