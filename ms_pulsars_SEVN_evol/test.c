#include <stdio.h>
#include <stdlib.h>

int main(){

	double temp[661380];
	FILE *temp_data=NULL;
	long np=0;
	temp_data=fopen("temp.txt","r");
	while(fscanf(temp_data,"%le\n",&temp[np])==1) {np++;}
	for(np=0;np<661380;np++){
		printf("%e\n",temp[np]);
	}
}
