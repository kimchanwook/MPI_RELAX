#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define COLS (5)
#define ROWS (16)

extern double get_change(double *x,double *y,int n);
extern void relax(double *dest,double *srce, int cols, int rows);
extern void init_grid(double **, double **, int cols, int rows);
extern void init_boundaries(double *,int, int);
extern void printf_buffer(double *, int, int);

int main(int argc, char const *argv[]){

  int size, rank;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Status status;
  int subrows=ROWS/size;
  double *local,*local_new, *p,*p_new;
  init_grid(&local,&local_new,COLS,subrows+2);

  if (0 == rank){
    init_grid(&p,&p_new,COLS,ROWS+2);
    init_boundaries(p,COLS,ROWS+2);
    memmove(p_new,p,COLS*(ROWS+2)*sizeof(double));
    printf_buffer(p_new,COLS,ROWS+2);
    printf("\n");}

  MPI_Barrier(MPI_COMM_WORLD);

  if(0 == rank)
    for(int j=0;j<size;j++) MPI_Send(p+j*subrows*COLS,(subrows+2)*COLS,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
  MPI_Recv(local, (subrows+2)*COLS, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  memmove(local_new,local,COLS*(subrows+2)*sizeof(double));
  MPI_Barrier(MPI_COMM_WORLD);

  int count=300;
  int up_number,dn_number;
  double *first_buff, *last_buff;
  init_grid(&first_buff,&last_buff,COLS,1);

  while(count-- >0){
    relax(local_new,local,COLS,subrows+2);
    memmove(last_buff,local_new,COLS*sizeof(double));
    memmove(first_buff,local_new+(subrows+1)*COLS,COLS*sizeof(double));
    up_number=(rank==size-1) ? MPI_PROC_NULL : rank+1;
    dn_number=(rank==0) ? MPI_PROC_NULL : rank-1;
    MPI_Sendrecv(local_new+(1)*COLS,COLS,MPI_DOUBLE,dn_number,0,first_buff,COLS,MPI_DOUBLE,up_number,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(local_new+(subrows)*COLS,COLS,MPI_DOUBLE,up_number,0,last_buff,COLS,MPI_DOUBLE,dn_number,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    MPI_Barrier(MPI_COMM_WORLD);
    memmove(local_new,last_buff,COLS*sizeof(double));
    memmove(local_new+(subrows+1)*COLS,first_buff,COLS*sizeof(double));
    memmove(local,local_new,COLS*(subrows+2)*sizeof(double));}

  MPI_Send(local_new+1*COLS,(subrows)*COLS,MPI_DOUBLE,0,(subrows+2),MPI_COMM_WORLD);

/*0 thread recieve data from local*/
  if (0==rank){
    for(int i=0;i<size;i++){
      MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      MPI_Recv(p_new+(status.MPI_SOURCE*subrows+1)*COLS,(subrows)*COLS,MPI_DOUBLE,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);}}

   MPI_Barrier(MPI_COMM_WORLD);

  if(0==rank){
    printf("\n");
    printf_buffer(p_new,COLS,ROWS+2);}

  free(local);
  free(local_new);
  free(p);
  free(p_new);
  free(first_buff);
  free(last_buff);
  MPI_Finalize();

return 0;}

void relax(double *pnew,double *pold, int cols, int rows){
  for ( int j = 1 ; j < (rows)-1; j++){
    for ( int i = 1 ; i < cols-1; i++){
      pnew[i+j*cols] = 0.25*(pold[i-1+j*cols]+pold[i+1+j*cols]+pold[i+(j-1)*cols]+pold[i+(j+1)*cols]);}}}

void init_boundaries(double *l_p,int cols,int rows){
  for (int i = 0 ; i < cols*rows; i++ )
  l_p[i] = 0;
  l_p[cols/2] = 1;}

void init_grid(double **p, double **p_new,int cols, int rows){
  if (NULL == (*p = malloc(cols * rows * sizeof(double)))){
    puts("Allocation Error.");
    exit(99);}

  if (NULL == (*p_new = malloc(cols*rows * sizeof(double)))) {
    puts("Allocation Error.");
    exit(99);}

  *p = (double *)malloc(cols*rows * sizeof(double));
  *p_new = (double *)malloc(cols*rows * sizeof(double));

  for ( int i = 0 ; i < cols*rows; i++){
    (*p)[i] = 0.;
    (*p_new)[i] = 0.;}}

void printf_buffer(double *p, int cols, int rows){
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      printf("%f\t",p[i*cols+j]);}
    puts(" ");}
  printf("\n");}

double get_change(double *x, double *y, int n){
  double result=0.;
  for (int i=0; i<n;i++)
    result+=x[i]*x[i]-y[i]*y[i];
return result;}

