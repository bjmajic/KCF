#include "jmyerror.h"
//错误处理重载函数：
void my_error_exit(j_common_ptr pInfo)
{
	my_error_mgrp myerr=(my_error_mgrp)pInfo->err;
  if(myerr && myerr->setjmp_buffer){
    longjmp(myerr->setjmp_buffer, 1);//longjmp执行到最后一个setjmp处。
  }
}
