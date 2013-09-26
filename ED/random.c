#include<stdio.h>
#include<string.h>
#include<stdlib.h>

void gen_random(char *s, const int len);

void main(int argc, char **argv){
	int len = atoi(argv[1]);
	char s1[len]; 
	char s2[len];
	FILE *out = fopen("input.txt", "w");
	gen_random(s1, len);
	gen_random(s2, len);
	fprintf(out, "%s %s", s1, s2);
} 
	
	


void gen_random(char *s, const int len) {
    static const char alphanum[] = "abc";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    s[len] = 0;
}

