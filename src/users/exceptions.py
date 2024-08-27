class CredentialsException(Exception):
    def __init__(self, message: str = None):
        if message:
            self.message = message
        else:
            self.message = "Could not validate credentials"


class NotEnoughPermissionException(Exception):
    def __init__(self, required):
        self.message = f"{required} permission is required"


class DuplicateEmailException(Exception):
    def __init__(self, email: str):
        self.email = email
        self.message = f"Email {email} already exists"


class EmailNotFoundException(Exception):
    def __init__(self, email: str):
        self.email = email
        self.message = f"Email {email} not found"
