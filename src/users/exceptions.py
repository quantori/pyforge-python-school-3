class CredentialsException(Exception):
    def __init__(self, message: str = None):
        if message:
            self.message = message
        else:
            self.message = "Could not validate credentials"
