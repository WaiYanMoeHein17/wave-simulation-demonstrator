"""
This is the exahype.toolkit package and holds the code of the "toolkit2", also called
"toolkit NG", the successor of the java-based ExaHyPE toolkit which itself is the
main glue code generator for ExaHyPE supposed to be runned by users.
"""

#Check Python version
from .configuration import checkPythonVersion
checkPythonVersion()

#High level API
from .configuration import checkDependencies
from .controller import Controller
