import copy

from pygromos.utils.typing import List, Dict, Iterable, Number

# FIELDS


class _generic_field:
    comment: str = ""

    fieldseperator = "\t"
    lineseperator = "\n"

    def __str__(self):
        return self.to_string()

    def __copy__(self):
        field = type(self)()
        for attr in vars(self):
            setattr(field, attr, getattr(self, attr))
        return field

    def to_string(self):
        raise NotImplementedError("to string method needs to be implemented!")


# BLOCKS
# genericblock:
class _generic_gromos_block:
    comment: str
    content: Iterable  # some content

    def __init__(self, name: str = None, used: bool = None, content: str = None):  # content:str,
        self.used = used
        self.name = name
        self.line_seperator = "\n"
        self.field_seperator = " \t "
        self.comment_char = "#"
        self._check_import_method(content=content)

    def __str__(self):
        return self.block_to_string()

    def __repr__(self):
        return str(self)

    def __iter__(self):
        return iter(self.content)

    def __copy__(self):
        block = type(self)(name=self.name, used=self.used, content=self.content)
        return block

    def __deepcopy__(self, memo):
        # return block as string, split by line and cut block title and END
        newContent = self.line_seperator.join(self.block_to_string().split(self.line_seperator)[1:-2])
        block = type(self)(content=newContent.split(self.line_seperator))
        return block

    def __eq__(self, __o: object) -> bool:
        if self.block_to_string() == str(__o):
            return True
        else:
            return False

    def _check_import_method(self, content: str = None):
        if content is not None:
            if isinstance(content, list) and all([isinstance(x, str) for x in content]):
                self.read_content_from_str(content)
            elif type(content) == self.__class__:
                self.content = content
            elif isinstance(content, str):
                self.read_content_from_str(content=content.split(self.line_seperator))
            else:
                raise Exception("Generic Block did not understand the type of content \n content: \n" + str(content))
        else:
            self.content = []

    def read_content_from_str(self, content: str):
        if isinstance(content, str):
            lines = content.split("\n")
        else:
            lines = content
        self.content = []
        for field in lines:
            if not field.strip().startswith("#") and not len(field.strip()) == 0:
                self.content.append(field.strip().split())

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        if (
            isinstance(self.content, list)
            and len(self.content) > 0
            and all([isinstance(x, (str, Number)) for x in self.content])
        ):
            result += self.field_seperator.join(map(str, self.content)) + self.line_seperator
        elif (
            isinstance(self.content, list)
            and len(self.content) > 0
            and all([isinstance(x, list) and all([isinstance(y, (str, Number)) for y in x]) for x in self.content])
        ):
            result += self.line_seperator.join(map(lambda x: self.field_seperator.join(map(str, x)), self.content))
        elif isinstance(self.content, (str, Number)):
            result += self.field_seperator + str(self.content) + self.line_seperator
        else:
            result += self.field_seperator + "EMPTY" + self.line_seperator
        result += self.line_seperator + "END" + self.line_seperator
        return result

    def get_name(self):
        return self.name


class _iterable_gromos_block(_generic_gromos_block):
    table_header = [""]

    def __init__(self, name: str, used: bool, content=None):
        self._content = []
        super().__init__(name, used, content=content)

    @property
    def content(self):
        return self._content

    @content.setter
    def content(self, content: Iterable):
        self._content = content

    def __getitem__(self, item: int):
        return self.content[item]

    def __len__(self):
        return len(self.content)

    def append(self, obj):
        self.content.append(obj)

    def extend(self, it: Iterable):
        self.content.extend(it)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += "#" + self.field_seperator + self.field_seperator.join(self.table_header) + "\n"
        for x in self.content:
            result += x.to_string()
        result += "END\n"

        return result


# GENERAL
class TIMESTEP(_generic_gromos_block):
    step: int
    t: float

    def __init__(
        self,
        t: float = None,
        step: int = None,
        content: Dict = None,
        subcontent: bool = False,
        name: str = "TIMESTEP",
        used: bool = True,
    ):

        if t is None and step is None:
            super().__init__(used=used, name=name, content=content)
        elif content is None:
            super().__init__(used=used, name=name, content=None)
            self.t = t
            self.step = step

        self.subcontent = subcontent

    def read_content_from_str(self, content: List[str]):
        content = content[0].strip().split()
        self.content = content
        self.t = float(content[1])
        self.step = float(content[0])

    def block_to_string(self) -> str:
        result = self.name + "\n"
        result += " \t{} \t{} \n".format(self.step, self.t)
        result += "END\n"
        return result


class TITLE(_generic_gromos_block):
    content: str
    field_seperator: str = "\t"
    line_seperator: str = "\n"
    order = [[["content"]]]
    pyGromosWatermark: str = ">>> Generated with PyGromosTools (riniker group) <<<"

    def __init__(
        self,
        content: str,
        field_seperator: str = "\t",
        line_seperator: str = "\n",
        name: str = "TITLE",
        used: bool = True,
    ):
        super().__init__(used=used, name=name, content=content)
        self.field_seperator = field_seperator
        self.line_seperator = line_seperator

    def __deepcopy__(self, memo):
        block = type(self)(content=None)
        block.content = copy.deepcopy(self.content)
        return block

    def read_content_from_str(self, content: List[str]):
        if type(content) == str:
            self.content = [content]
        else:
            self.content = content

    def block_to_string(self) -> str:
        result = ""
        result += str(self.name) + self.line_seperator
        result += "".join(self.content)
        if self.pyGromosWatermark not in result:
            result += self.line_seperator + self.field_seperator + self.pyGromosWatermark + self.line_seperator
        result += "END" + self.line_seperator
        return result


class TRAJ(_iterable_gromos_block):
    content: Iterable

    def __init__(self, timestep_blocks: Iterable):
        super().__init__(used=True, name="Trajectory")
        self.content = timestep_blocks
        self.dt = 0

    def block_to_string(self) -> str:
        return self.name + " contains \t" + str(len(self.content))
