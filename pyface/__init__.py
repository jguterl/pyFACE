print('import face ---')import sys, ossys.setdlopenflags(os.RTLD_LAZY | os.RTLD_GLOBAL)from . import *#from . import FACECfrom .FACEC import *#from .pyface import *from Forthon import *#from . import facepy# Order in which packages are called matters due to dependencies (added by J.Guterl)from .facepy import *from . import *__version__ = '1.0.1'from .__git__ import sha, branch ,tagprint('>>> FACE version:',__version__) print('>>> FACE git:[{}][{}][{}]'.format(sha[:6],branch, tag))    