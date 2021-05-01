#include <sys/socket.h>
#include <sys/un.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "main.h"


char *socket_path = "/tmp/socket";

#define BUF_SIZE 4096			/* block transfer size */

int main(int argc, char **argv)
{
	int c, uds, bytes,ngenomas,NCHILDS,i,j;
	long time_usec_begin;
	long time_usec_end;
	long elapsed_time;
    char buf[BUF_SIZE];			/* buffer for incoming file */
    struct sockaddr_un channel;		/* Unix Domain socket */

	if (argc != 3) {
		printf("Usage: client file-name\n");
		exit(1);
	}

	if ( (uds = socket(AF_UNIX, SOCK_STREAM, 0)) == -1) {
		perror("socket error");
		exit(1);
	}

	memset(&channel, 0, sizeof(channel));
	channel.sun_family= AF_UNIX;
	strncpy(channel.sun_path, socket_path, sizeof(channel.sun_path)-1);

	if (connect(uds, (struct sockaddr*)&channel, sizeof(channel)) == -1) {
		perror("connect error");
		exit(1);
	}
    /* Connection is now established. Send file name including 0 byte at end. */

	NCHILDS=atoi(argv[1]);
	i=atoi(argv[2]);

	get_time_useconds(&time_usec_begin);


	GENOME_LIST * gl = (GENOME_LIST *) calloc(1, sizeof(GENOME_LIST));
	read_genomes(gl, "input/cds.fna");
	printf("Filho %d Read: %ld genomes\n",i, gl->n_genomes);
	ngenomas=gl->n_genomes;


	MUTATION_ARRAY * mutation_array = (MUTATION_ARRAY *) calloc(1, sizeof(MUTATION_ARRAY));
	GENOME * g = gl->phead;

	int x;
	x=getpid();

	for (j = 0; j < i; j++)
	{

		g=g->pnext;
	}


	while (g != NULL) {
		genome_cmp(g, mutation_array);
		send_mutation_array(mutation_array,uds, x, 0);
		free_mutations(mutation_array);
		for (int k = 0; k < NCHILDS; k++)
		{
			g=g->pnext;
		}
		kill(getppid(),SIGUSR1);
	}

	close(uds);
	exit(0);


}


void send_mutation_array (MUTATION_ARRAY * mutation_array,int uds, int getpid, int detail) {
	char * buf = (char *) malloc(sizeof(char) * 1000000);
	for (int i=0 ; i<mutation_array->n_mutations ; i++) {
		MUTATION * aux = mutation_array->mutations + i;
		sprintf(buf, "#%d|%s;%s|%s|(%d)\n", getpid ,aux->genome_a, aux->genome_b, aux->gene, aux->seq_mutations.n);
		if (detail) {
			for (int j = 0; j < aux->seq_mutations.n; j++) {
				sprintf(buf, "%s%d;", buf, aux->seq_mutations.arr[j]);
			}
		}

		writen(uds,buf,strlen(buf));
	}
}

/**
	 * Frees a given mutation array
	 * @param mutation_array - pointer to the previously allocated mutation array
	 */
	void free_mutations(MUTATION_ARRAY * mutation_array) {
		for (int i=0 ; i<mutation_array->n_mutations ; i++) {
			free((mutation_array->mutations + i)->seq_mutations.arr);
		}

		free(mutation_array->mutations);

		mutation_array->n_mutations = mutation_array->size_mutations = 0;
		mutation_array->mutations = NULL;
	}

ssize_t readn(int fd, void *vptr, size_t n) {
	size_t  nleft;
	ssize_t nread;
	char   *ptr;

	ptr = vptr;
	nleft = n;
	while (nleft > 0) {
		if ( (nread = read(fd, ptr, nleft)) < 0) {
			if (errno == EINTR)
           nread = 0;      /* and call read() again */
				else
					return (-1);
			} else if (nread == 0)
             break;              /* EOF */

			nleft -= nread;
			ptr += nread;
		}
     return (n - nleft);         /* return >= 0 */
	}

ssize_t writen(int fd, const void *vptr, size_t n)
	{
		size_t nleft;
		ssize_t nwritten;
		const char *ptr;

		ptr = vptr;
		nleft = n;
		while (nleft > 0) {
			if ( (nwritten = write(fd, ptr, nleft)) <= 0) {
				if (nwritten < 0 && errno == EINTR)
                 nwritten = 0;   /* and call write() again */
					else
                 return (-1);    /* error */
				}
			nleft -= nwritten;
			ptr += nwritten;
		}
		return (n);
	}


	/**
	 * Inserts a new genome at the tail of given genome list
	 * @param gl - genome list
	 * @param g - new genome
	 */
void insert_genome (GENOME_LIST * gl, GENOME * g) {
	g->pnext = NULL;
	g->pprev = NULL;

	if (gl->phead == NULL)
		gl->phead = g;

	if (gl->ptail != NULL) {
		g->pprev = gl->ptail;
		gl->ptail->pnext = g;
	}

	gl->ptail = g;
	gl->n_genomes++;
}

	/**
	 * Searches for a gene in a given genome
	 * @param genome - genome to be scanned
	 * @param gene_name - gene to searched
	 * @return - pointer to the found gene or NULL if no match
	 */
GENE * find_gene (GENOME * genome, char * gene_name) {
	for (int i=0 ; i<genome->n_genes ; i++) {
		if (strcmp((genome->genes+i)->name, gene_name) == 0)
			return genome->genes+i;
	}
	return NULL;
}

	/**
	 * Inserts a new element into a given int_array
	 * @param int_array - given integer array
	 * @param element - element to be inserted
	 */
void insert_int_array (INT_ARRAY * int_array, int element) {
	if (int_array->n >= int_array->size) {
		int_array->size = (int_array->size != 0) ? int_array->size * 2 : 2;
		int_array->arr = (int *) realloc(int_array->arr, int_array->size * sizeof(int));
	}

	int_array->arr[int_array->n] = element;
	int_array->n++;
}

	/**
	 * Compares two genes by subtracting each of the nucleotide sequences values of g1 and g2
	 * @param g1 - gene 1 to be compared
	 * @param g2 - gene 2 to be compared
	 * @return integer array containing the differences between the tow genes
	 */
INT_ARRAY * gene_cmp (GENE g1, GENE g2) {
	int i;
	INT_ARRAY * to_return = (INT_ARRAY *) calloc(1, sizeof(INT_ARRAY));

	for (i=0 ; *(g1.seq+i) != '\0' ; i++) {
		int x = abs((int) *(g1.seq+i) - (int) *(g2.seq+i));
		if (x != 0) insert_int_array(to_return, i);
	}

	return to_return;
}

	/**
	 * Inserts a new mutation into a given mutation array
	 * @param mutation_array - array of mutations
	 * @param genome_a - genome used for comparison against genome b
	 * @param genome_b - genome used for comparison against genome a
	 * @param gene - gene on which the two genomes were previously compared
	 * @param gene_mut - integer array with all the found mutations
	 */
void insert_mutation (MUTATION_ARRAY * mutation_array, char * genome_a, char * genome_b, char * gene, INT_ARRAY * gene_mut) {
	if (mutation_array->n_mutations >= mutation_array->size_mutations) {
		mutation_array->size_mutations = (mutation_array->size_mutations != 0) ? mutation_array->size_mutations * 2 : 2;
		mutation_array->mutations = (MUTATION *) realloc(mutation_array->mutations, mutation_array->size_mutations * sizeof(MUTATION));
	}
	MUTATION * aux = mutation_array->mutations + mutation_array->n_mutations;
	strcpy(aux->genome_a, genome_a);
	strcpy(aux->genome_b, genome_b);
	strcpy(aux->gene, gene);

	aux->seq_mutations = *gene_mut;

	mutation_array->n_mutations++;
}

	/**
	 * Compares a given genome against all its subsequent genemoes in a genome list
	 * @param genome - reference genome to compare against all the subsequent genomes
	 * @param mutation_array - array in which the comparison results (mutations) will be stored
	 */
void genome_cmp (GENOME * genome, MUTATION_ARRAY * mutation_array) {
	GENE * base_gene;
	INT_ARRAY * gene_mut = NULL;

	for (int i=0 ; i<genome->n_genes ; i++) {
		base_gene = genome->genes+i;

		GENOME * tmp_genome = genome->pnext;
		while (tmp_genome != NULL) {
			GENE * new_gene = find_gene (tmp_genome, base_gene->name);
			if (new_gene != NULL) {
				if ((gene_mut = gene_cmp(*base_gene, *new_gene)) != NULL) {
					insert_mutation(mutation_array, genome->name, tmp_genome->name, base_gene->name, gene_mut);
				}
			}
			tmp_genome = tmp_genome->pnext;
		}
	}
}

	/**
	 * Removes white spaces ' ' and '\n' from a given sting
	 * @param str - string with no ' ' or '\n'
	 */
void remove_white_spaces (char * str) {
	int c=0, j=0;
	while(str[c] != '\0') {
		if(str[c] != ' ' && str[c] != '\n')
			str[j++] = str[c];
		c++;
	}
	str[j]='\0';
}

	/**
	 * Searches, by name, for a given gene in a known gene dictionary
	 * @param name - gene name to search for
	 * @return - pointer to the found dictionary entry or NULL if non-existent
	 */
GENE_DICT * find_gene_dict (char * name) {
	for (int i=0 ; i<DICT_SIZE ; i++)
		if (strcmp(name, gd[i].name) == 0)
			return gd+i;
		return NULL;
	}

	/**
	 * Finds the number of dummy nucleotides to append to the nucleotide sequence
	 * Not required but useful if displaying a square matrix with the gene comparison result
	 * @param name - gene name
	 * @return - number of dummy nucleotides to append
	 */
	int get_gene_padding (char * name) {
		GENE_DICT * gene = find_gene_dict (name);
		if (gene != NULL)
			return gene->padding;
		return 0;
	}

	/**
	 * Creates a new gene given a gene name and a nucleotide sequence
	 * @param name - new gene name
	 * @param seq - new gene nucleotide sequence
	 * @return - pointer to the created gene
	 */
	GENE * create_gene (char * name, char * seq) {
		GENE * ret = (GENE *) malloc(sizeof(GENE));
		remove_white_spaces(seq);
		int N = get_gene_padding(name);
		ret->seq = (char *) malloc(sizeof(char) * (strlen(seq) + N + 1));

		strcpy(ret->name, name);
		sprintf(ret->seq, "%s%.*s", seq, N, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
		return ret;
	}

	/**
	 * Inserts a new gene into a given genome
	 * @param genome - pre-existing genome
	 * @param gene - gene to be inserted in the given genome
	 */
	void insert_gene (GENOME * genome, GENE * gene) {
		if (genome->n_genes >= genome->size_genes) {
			genome->size_genes = (genome->size_genes != 0) ? genome->size_genes * 2 : 2;
			genome->genes = (GENE *) realloc(genome->genes, genome->size_genes * sizeof(GENE));
		}

		GENE * g = genome->genes + genome->n_genes;
		*g = *gene;

		genome->n_genes++;
	}

	void read_genomes (GENOME_LIST * gl, char * path) {
		long bytes, total=0, size;

		char * cds=NULL;

		int fd = open(path, O_RDONLY);
		if (fd == -1) { perror("File open"); exit(1); }

		size = lseek(fd, 0, SEEK_END);
		lseek(fd, 0, SEEK_SET);

		cds = (char *) malloc(sizeof(char) * (size+1));
		while ((bytes = read(fd, cds+total, BUF_SIZE)))
			total += bytes;

		close(fd);

		parse_genome(gl, cds);
	}

	char * find_protein_name (char * protein) {
		for (int i=0 ; i<DICT_SIZE ; i++) {
			if (strcmp(protein, gd[i].prot) == 0)
				return gd[i].name;
		}
		return "";
	}

	GENOME * find_genome (GENOME_LIST * gl, char * g_id) {
		if (gl == NULL) return NULL;
		if (gl->phead == NULL || gl->ptail == NULL) return NULL;

		GENOME * to_return = gl->ptail;
		while (to_return != NULL) {
			if (strcmp(g_id, to_return->name) == 0)
				return to_return;
			to_return = to_return->pprev;
		}

		return NULL;
	}

	/**
	 * Parses a given code region sequence by genomes and genes,
	 * populating the received genome list with the loaded values
	 * @param gl - pointer to the genome list
	 * @param cds - loaded given code region sequence containing all the genomes
	 */
	void parse_genome (GENOME_LIST * gl, char * cds) {
		int n=0;

		char * token;
		char needle[] = ">";
		char genome_id[MAX100], protein[MAX100];

		token = strtok(cds, needle);

		while (token != NULL) {
			sscanf(token, "%[^.]%*s%s%*[^\n]%n", genome_id, protein, &n);

			strcpy(protein, find_protein_name(protein));
			if (strcmp(protein, "") != 0) {

				GENE *new_gene = create_gene(protein, token + n + 1);

				GENOME * p_genome = find_genome(gl, genome_id);
				if (p_genome == NULL) {
					p_genome = (GENOME *) calloc(1, sizeof(GENOME));
					strcpy(p_genome->name, genome_id);
					insert_genome(gl, p_genome);
				}
				insert_gene(p_genome, new_gene);
			}

			token = strtok(NULL, needle);
		}

		free(cds);
	}

	/**
	 * prints a given genome to the std output
	 * @param genome - genome to be printed
	 */
	void print_genome (GENOME genome) {
		GENE * gene = genome.genes;

		printf("Genome: %s, %d\n", genome.name, genome.n_genes);
		for (int i=0 ; i<genome.n_genes ; i++) {
			printf("\tName: %s\n", gene->name);
			printf("\tSequence: %s\n\n", gene->seq);

			gene++;
		}
	}

	/**
	 * Gets the number of microseconds elapsed since 1/1/1970
	 * @param time_usec - variable in which the elapsed time will be stored
	 * @return - elapsed time since 1/1/1970
	 */
	long get_time_useconds(long * time_usec) {
		struct timeval time;
		gettimeofday(&time, NULL);

		*time_usec = (long) (time.tv_sec * 1000000 + time.tv_usec);
		return *time_usec;
	}
