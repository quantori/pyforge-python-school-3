from sqlalchemy.orm import Mapped
from database import Base, str_uniq, int_pk


class Molecule(Base):
    id: Mapped[int_pk]
    name: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.id}, "
            f"name={self.name!r},"
        )

    def __repr__(self):
        return str(self)
