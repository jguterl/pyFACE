/* Created by Ron Cohen */
/* $Id: uedgeC_Forthon.c,v 7.1 2019/11/01 22:20:49 meyer8 Exp $ */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of UEDGE. This is needed since the modules      */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include <Python.h>
#define NPY_NO_DEPRECATED_API 8
#include <numpy/arrayobject.h>


#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

static PyObject *ErrorObject;


/* ######################################################################### */
/* # Method list                                                             */
static struct PyMethodDef FACEC_methods[] = {
  {NULL,NULL}};


static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "FACEC", /* m_name */
  "FACEC", /* m_doc */
  -1,                  /* m_size */
  FACEC_methods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
  };




/* ######################################################################### */
/* # The initialization function                                             */
PyMODINIT_FUNC PyInit_FACEC(void)

{

  PyObject *m, *d;
  m = PyModule_Create(&moduledef);


  d = PyModule_GetDict(m);
  ErrorObject = PyErr_NewException("FACEC.error",NULL,NULL);
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module FACEC");

  /*pystdout = PySys_GetObject("stdout");*/
  /*PyFile_WriteString("Forthon edition\n",pystdout);*/

  import_array();


  return m;


}


