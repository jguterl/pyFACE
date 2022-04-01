/**
 * @file exmain.c
 *
 * Purpose: Allow the trapping of control-C to support Basis-like debug
 *          mode in the python version.
 *
 * $Id: exmain.c,v 1.4 2019/12/19 16:44:25 meyer8 Exp $
 *
 */


#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#ifdef HAS_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif



static struct sigaction act,oact;
static sigjmp_buf ev;


/*
    Handler for SIGINT signal
*/
void int_handler() {
   char mymyline[200],*ret;
   sigset_t block_mask;
   printf("\nType \"cont\" to continue exmain(), \"abort\" (not compatible with openmp) or \"stop\" (with openmp) to return to Python prompt \n");
   printf("or a single line to be evaluated by Python.\n");
#pragma omp master
    {
int condition;
condition=1;
   while(1){


       printf("Debug>>> ");
       ret = fgets(mymyline,150,stdin);
       if(ret == (char *)NULL)return;

       if(strncmp(mymyline,"cont",4) == 0){
           return;
       } else if (strncmp(mymyline,"abort",5) == 0) {
          PyRun_SimpleString("bbb.exmain_aborted = True");
          siglongjmp(ev,1);
       } else if (strncmp(mymyline,"stop",4) == 0) {
	  PyRun_SimpleString("print(\"Stopping exmain ... Please wait...\")");
          PyRun_SimpleString("bbb.exmain_aborted = True");
          return;
       } else if (strncmp(mymyline,"exit",4) == 0) {
          PyRun_SimpleString("bbb.exmain_aborted = True");
          siglongjmp(ev,1);
       } else {
          PyRun_SimpleString(mymyline);
          /* matplotlib seems to unset the hander
             so it is set again just in case */
          sigfillset(&block_mask);
          act.sa_handler = int_handler;
          act.sa_mask = block_mask;
          act.sa_flags = 0;
          sigaction(SIGINT,&act,NULL);
       }
   }

}
}




/* FORTHON is defined by the Python build. This exmain does nothing when
   compiled for the basis version of the code, it just drops through to
   the Fortran routine.  */


void initializec_() {


   sigset_t block_mask;
   int ival;
#pragma omp master
{
   ival = sigsetjmp(ev,1);
   if(ival != 0){
       sigaction(SIGINT,&oact,NULL);
       return;
   }
   }


/* setup to catch SIGINT and save the previous handler to be restored
   on return */
#pragma omp master
{
   sigfillset(&block_mask);
   act.sa_handler = int_handler;
   act.sa_mask = block_mask;
   act.sa_flags = 0;
   sigaction(SIGINT,&act,&oact);

   }





/*  now call the Fortran version of exmain  */

 


   sigaction(SIGINT,&oact,NULL);

}




