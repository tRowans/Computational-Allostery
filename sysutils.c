#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>
#include <ctype.h>
#include "lib.h"

//External program calling

int spawn(char* program,char** arg_list)
//runs program with arguments arg_list and returns processid
{
    pid_t child_pid;
    int child_status;
    child_pid=fork();//duplicates processs
    if(child_pid!=0)//if parent
    {
        wait(&child_status);//pauses parent until child returns
        if(WIFEXITED(child_status))
        {
            printf("\nchild exited normally with exit code %d\n\n",WEXITSTATUS(child_status));
        }
        else
        {
            printf("\nchild exited abnormally\n");
            exit(0);
        }
        return child_pid;
    }
    else//executes program in child
    {
        execvp(program,arg_list);
        fprintf(stderr,"error occured in execvp with %s\n",program);
        abort();
    }
}

void copy(char path[],char dpath[])
//copy file path to dpath
{
    char* arg_list[]={"cp","-f",path,dpath,NULL};
    spawn("cp",arg_list);
}

void mv(char path[],char dpath[])
//move file path to dpath
{
    char* arg_list[]={"mv","-f","-T",path,dpath,NULL};
    spawn("mv",arg_list);
}

void ddpt(char fname[])
//call DDPT
{
    printf("\n\nCalling DDPT\n\n");

    char* arg_list1[]={"genENM","-pdb",fname,"-het","-ca","-c","8",NULL};//list of arguments for spawn
    spawn("/usr/local/DDPT/GENENMM",arg_list1);

    char* arg_list2[]={"diagstd",NULL};//list of arguments for spawn
    spawn("/usr/local/DDPT/DIAGSTD",arg_list2);

    char* arg_list3[]={"freqeuncy",NULL};//list of arguments for spawn
    spawn("/usr/local/DDPT/FREQEN",arg_list3);
}

void enms(int d)
{
    char cpy[50];
    sprintf(cpy,"runs/%d.vmd",d);
    copy("ENM.vmd",cpy);
    printf("\n%d ENM file stored for animation\n",d);
}

void cpdir(char src[],char dest[])
{
	char* arg_list[]={"cp","-f","-r",src,dest,NULL};
    spawn("cp",arg_list);
}

void rmvdir(char src[])
{
	char* arg_l[]={"rm","-r","-f",src,"NULL"};
	spawn("rm",arg_l);
}
//----------------------------
