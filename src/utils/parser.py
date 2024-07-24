from fastapi import HTTPException, UploadFile


async def parse_text_file(file: UploadFile):
    if file.content_type != "text/plain":
        raise HTTPException(status_code=422, detail="Incorrect file format")

    content = await file.read()
    content_str = content.decode("utf-8")
    cleaned_smiles = [
        smile.strip().replace("\n", "").replace("\r", "")
        for smile in content_str.split(",")
    ]

    return cleaned_smiles


if __name__ == "__main__":
    pass
