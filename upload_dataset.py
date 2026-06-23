from huggingface_hub import HfApi

REPO_ID = "FAKIrik/Seismo_datasets"
REPO_TYPE = "dataset"
LOCAL_FOLDER = "./1st_selection"

api = HfApi()

# Создаём репозиторий, если его нет
api.create_repo(repo_id=REPO_ID, repo_type=REPO_TYPE, exist_ok=True)

# Загружаем всю папку в корень репозитория
api.upload_large_folder(
    repo_id=REPO_ID,
    repo_type=REPO_TYPE,
    folder_path=LOCAL_FOLDER,
)

print("✅ Загрузка завершена!")
