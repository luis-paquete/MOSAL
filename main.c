#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#define GAPS 	 "/Gaps"
#define INDELS 	 "/Indels"

#define DP 		 "/DP"
#define DPP		 "/DP_prune"

#define SS		 "/Subs_score"
#define MAT		 "/Matches"

#define NO_TB	 "/No_traceback"
#define TB		 "/Traceback"

void choose_program( char prog[], char *args[], int argc, char const * argv[] );
int check_main_args( int argc, char const *argv[] );

int main(int argc, char const *argv[]){

	char prog[512]="";		// paths to external files can be big
	char *args[6];


	printf("%s\n", argv[1]);
	
	if( argc < 5 || !check_main_args(argc, argv) ){
		printf("\nUsage: %s seq1_file seq2_file [gaps|indels] [dp|dpp -b=NUMBER] [-ss=FILE] [--no-traceback]\n", argv[0]);
		printf("\nOptions:\n");

		printf("  seq1_file\t\tpath to the 1st sequence file (FASTA format)\n");
		printf("  seq2_file\t\tpath to the 2nd sequence file (FASTA format)\n");

		printf("  gaps\t\t\tspecify the problem to be relative to gaps.\n");
		printf("  indels\t\tspecify the problem to be relative to indels.\n");
		
		printf("  dp\t\t\tuse Dynamic Programming approach.\n");
		printf("  dpp\t\t\tuse Dynamic Programming approach with pruning (need to specify -b=NUMBER).\n");

		printf("  -ss=FILE\t\tuse substitution score(ss) instead of #matches (path to the ss table).\n");
		printf("  -b=NUMBER\t\tspecify the number of MID bounds to prune.\n");
		printf("  --no-traceback\toutput only the scores without the result alignments.\n");

		return -1;
	}

	// relative path to program from the location execution started
	char path_to_prog[300];
	char * last_slash = strrchr(argv[0], '/');
	memcpy( path_to_prog, argv[0], last_slash-argv[0] ) ;

	// add relative path to current program
	strcat(prog, path_to_prog);

	choose_program( prog, args, argc, argv );
	printf("running \"%s ", prog);

	int i;
	for( i = 1; i < 5; ++i )
		if(args[i])
			printf("%s ", args[i]);
	printf("\"\n");

	if (execv(prog, args) == -1)
    	perror("Error executing program");


	return 0;
}

int check_main_args( int argc, char const *argv[] ){
	return (  ( strcmp(argv[3],"gaps")==0 || strcmp(argv[3],"indels")==0 ) && 
			  ( strcmp(argv[4],"dp")==0 || strcmp(argv[4],"dpp")==0 )
		   );
}

void choose_program( char prog[], char *args[], int argc, char const * argv[] ){
	int i;
	char * ss = NULL;
	int tb = 1;
	int num_args = 3;

	args[0] = "prog";
	args[1] = (char *) argv[1];
	args[2] = (char *) argv[2];
	args[3] = NULL;
	args[4] = NULL;
	args[5] = NULL;


	// load sequences




	// gaps or indels
	if( !strcmp( argv[3], "gaps" ) )
		strcat( prog, GAPS );
	else
		strcat( prog, INDELS );
	


	// match or subs_score
	for( i = 3; i < argc && !ss ; ++i )
		ss=strstr( argv[i], "-ss=" );

	if(ss){
		strcat( prog, SS );
		args[num_args++] = strchr(ss, '=')+sizeof(char);	// keep just the value of '-ss'
	}
	else
		strcat( prog, MAT );


	// traceback or no traceback
	for( i = 3; i < argc && tb; ++i )
		tb=strcmp(argv[i], "--no-traceback");
	
	if(!tb) strcat( prog, NO_TB);
	else	strcat( prog, TB);


	// choose algorithm
	if( !strcmp( argv[4], "dp" ) )
		strcat(prog, DP);
	else if( !strcmp( argv[4], "dpp" ) ){
		strcat(prog, DPP);
		// number of bounds
		char * b = NULL;
		for( i = 3; i < argc && !b; ++i )
			b=strstr( argv[i], "-b=" );
		
		if(b)
			args[num_args] = strchr(b,'=')+sizeof(char);	// keep just the value of '-b'
		else{
			printf("Prune version specified but no number of bounds (-b=X)\n");
			exit(-1);
		}

	}
	else{
		printf("Wrong algorithm chosen: %s\n", argv[4]);
		exit(-1);
	}

	strcat(prog, "/prog");

}

