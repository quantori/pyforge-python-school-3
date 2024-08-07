import argparse
from fastapi import FastAPI
from pathlib import Path

app = FastAPI()

# Global variable to store the mount directory
mount_dir = None


@app.get("/")
def read_root():
    return {"message": f"Mounted at {mount_dir}"}


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run a FastAPI application with a required --mount-dir argument.")
    parser.add_argument("--mount-dir", type=str, required=True, help="The directory to mount the database files. "
                                                                     ".collection directory will be created under here.")

    # Parse arguments
    args = parser.parse_args()

    # Validate the mount directory
    global mount_dir
    mount_dir = args.mount_dir
    if not Path(mount_dir).exists():
        print(f"Error: The specified mount directory '{mount_dir}' does not exist.")
        exit(1)

    # Run the FastAPI app
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)





if __name__ == "__main__":
    main()
