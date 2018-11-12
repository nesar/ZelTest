#include<stdio.h>
/* Our structure */
struct record
{ 
    int a,b,c;
};


int main()
{
  int count;
  FILE *ptr;
  struct record myRecord;
  ptr=fopen("test.bin","rb");
  if (!ptr)
      { 
        printf("Unable to open file!");     return 1;
      }
  for ( count=1; count <= 10; count++)
      { 
        fread(&myRecord,sizeof(struct record),1,ptr);   printf("%d\n",myRecord.a);
      }  fclose(ptr);
      return 0;
}
