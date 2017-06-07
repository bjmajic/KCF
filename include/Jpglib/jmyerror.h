#include "jinclude.h"
#include "jpeglib.h"
#include "jerror.h"
#include <setjmp.h>

typedef struct my_error_mgr{
  struct jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
}my_error_mgr, *my_error_mgrp;

void my_error_exit(j_common_ptr pInfo);
