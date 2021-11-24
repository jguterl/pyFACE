#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys
import os
import os.path
import string
import site
from Forthon.compilers import FCompiler
import getopt

GitHash=''
GitRemoteRepo=''
GitBranch=''
GitTag=''
UEDGEfolder=os.getcwd()
GitRepo=''
try:
    import git #gitpython
    Repo=git.Repo()
    GitHash=Repo.head.object.hexsha
    GitBranch=Repo.active_branch.name
    #GitRepo=Repo.active_branch.repo.name
except:
    pass
# Getting version from git tag
#try:
    #version=Repo.tags[-1].name
#except:
version='1.0'


try:
    os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'
    import distutils
    from distutils.core import setup
    from distutils.core import Extension
    from distutils.dist import Distribution
    from distutils.command.build import build
    from subprocess import call
    import numpy
except:
    raise SystemExit("Distutils problem")

optlist, args = getopt.getopt(sys.argv[1:], 'gt:F:', ['parallel', 'petsc','omp','nersc'])
machine = sys.platform
debug = 0
fcomp = None
parallel = 0
petsc = 0
OMP=False
for o in optlist:
    if o[0] == '-g':
        debug = 1
    elif o[0] == '-t':
        machine = o[1]
    elif o[0] == '-F':
        fcomp = o[1]
    elif o[0] == '--parallel':
        parallel = 1
    elif o[0] == '--petsc':
        petsc = 1
    elif o[0] == '--omp':
        OMP = True
    elif o[0] == '--fcomp':
        fcomp=o[1]

# OMP add-on
#OMPpackages=['bbb','com','api']
#OMPlisthtreadprivatevars='../../ppp/ListVariableThreadPrivate_final.txt'
CARGS=[]
FARGS=['-g -fmax-errors=15', '-DFORTHON','-cpp','-Wconversion','-fimplicit-none']
if OMP:
    FARGS=FARGS+['-fopenmp']
    CARGS=CARGS+['-fopenmp']
    OMPargs=['--omp']
else:
    OMPargs=[]
OMPFLAGS='OMPFLAGS = {}'.format(' '.join(OMPargs))

# Flags for makefile. Flags are easier to handle from setup.py and it prevents dealing with the makefile.)

FARGSDEBUG=['-fbacktrace','-ffree-line-length-0', '-fcheck=all','-ffpe-trap=invalid,overflow,underflow -finit-real=snan','-Og']
FARGSOPT=['-Ofast']

if debug==1:
    FARGS=FARGS+FARGSDEBUG
else:
    FARGS=FARGS+FARGSOPT

FLAGS ='DEBUG = -v --fargs "{}"'.format(' '.join(FARGS))
if CARGS!=[]:
    FLAGS =FLAGS+' --cargs="{}"'.format(' '.join(CARGS))



sys.argv = ['setup2.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)


class uedgeBuild(build):
    def run(self):

        call(['make',FLAGS,OMPFLAGS, '-f', 'Makefile.Forthon3'])
        build.run(self)


class uedgeClean(build):
    def run(self):
        call(['make', '-f', 'Makefile.Forthon3', 'clean'])


facepkgs = ['face']


def makeobjects(pkg):
    return [pkg+'_p.o', pkg+'pymodule.o']


uedgeobjects = []

# add here any extra dot o files other than pkg.o, pkg_p.o


dummydist = Distribution()
dummydist.parse_command_line()
dummybuild = dummydist.get_command_obj('build')
dummybuild.finalize_options()
builddir = dummybuild.build_temp

uedgeobjects = map(lambda p: os.path.join(builddir, p), uedgeobjects)
library_dirs = fcompiler.libdirs
libraries = fcompiler.libs


with open('pyscripts/__version__.py','w') as ff:
    ff.write("__version__ = '%s'\n"%version)
    ff.write("GitTag='{}'\n".format(GitTag))
    ff.write("GitBranch='{}'\n".format(GitBranch))
    ff.write("GitHash='{}'\n".format(GitHash))

define_macros=[("WITH_NUMERIC", "0"),
               ("FORTHON_PKGNAME", '\"uedgeC\"'),
               ("FORTHON","1")]

# check for readline
rlncom = "echo \"int main(){}\" | gcc -x c -lreadline - "
rln = os.system(rlncom)
if rln == 0:
   define_macros = define_macros + [("HAS_READLINE","1")]
   os.environ["READLINE"] = "-l readline"
   libraries = ['readline'] + libraries


setup(name="face",
      version=version,
      author='R. Smirnov/J. Guterl',
      author_email="",
      maintainer='',
      maintainer_email='',
      description="FACE",
      platforms="Unix, Windows (cygwin), Mac OSX",
      packages=['face'],
      package_dir={'face': 'pyscripts'},
      # include_package_data=True,
      scripts=[],
      ext_modules=[Extension('face.FACEC',
                             ['FACEC_Forthon.c',
                              os.path.join(builddir, 'Forthon.c'),'face/handler.c','face/initialize.c'],
                             include_dirs=[builddir, numpy.get_include()],
                             library_dirs=library_dirs,
                             libraries=libraries,
                             define_macros=define_macros,
                             extra_objects=uedgeobjects,
                             extra_link_args=CARGS+['-g','-DFORTHON'] +
                             fcompiler.extra_link_args,
                             extra_compile_args=fcompiler.extra_compile_args
                             )],

      cmdclass={'build': uedgeBuild, 'clean': uedgeClean},
      install_requires=['forthon'],
      # note that include_dirs may have to be expanded in the line above
      classifiers=['Programming Language :: Python :: 3']
      )
