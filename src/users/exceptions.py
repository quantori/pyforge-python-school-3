class CredentialsException(Exception):
    def __init__(self, message: str = None):
        if message:
            self.message = message
        else:
            self.message = "Could not validate credentials"


class DuplicateEmailException(Exception):
    def __init__(self, email: str):
        self.email = email
        super().__init__(f"Email {email} is already registered")


class EmailNotFoundException(Exception):
    def __init__(self, email: str):
        self.email = email
        super().__init__(f"Email {email} not found")
