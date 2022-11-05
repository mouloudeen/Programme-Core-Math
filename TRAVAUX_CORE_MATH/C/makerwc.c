#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>

int main(){
FILE *fic1 = fopen("log1.wc", "r");
if (fic1 == NULL)
    exit(1);

FILE *fic2 = fopen("log2.wc","r");
if (fic2 == NULL)
    exit(1);

FILE *fic3 = fopen("log.wc","w");
if (fic3 == NULL)
    exit(1);
double x;

while (!feof(fic1)){
        
        fscanf(fic1,"%la",&x);
        fprintf(fic3,"%la\n",x);
}
printf("log1.wc finis\n");
while (!feof(fic2)){
        
        fscanf(fic2,"%la",&x);
        fprintf(fic3,"%la\n",x);
}
fclose(fic1);
fclose(fic2);
fclose(fic3);
return 0;
}
/*gcc -o makerwc.o makerwc.c
./makerwc.o*/