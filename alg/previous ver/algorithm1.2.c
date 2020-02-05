#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mkl.h>


int main(int argc, char *argv[] ) {
  FILE *fp;
  int *row_ptr,*fill,*col_ind;
  double *val,a;
  char tmp[256];
  int i,j;
  int *nnonzero_row;
  int n,nnonzero;
  int row,col;

  //Usage error
  if(argc !=2){
    printf("Usage: sample <input_filename>\n");
    exit(1);
  }
  //file open error
  if((fp=fopen(argv[1],"r"))==NULL){
    printf("file open error!\n");
    exit(1);
  }

  i=0;//nnonzero counter
  //read file
  while( fgets(tmp,sizeof(tmp),fp) ){
    if(tmp[0] == '%'){
      //ignore commments
    }else{
      if(i == 0){
        sscanf(tmp,"%d %d %d", &n, &n, &nnonzero);//number of row,col,nonzero
        // printf("n:%d n:%d nnonzero:%d\n",n,n,nnonzero);
        nnonzero_row = (int *)malloc(sizeof(int)*n);
        for(j=0;j<n;j++){//initialize
          nnonzero_row[j] = 0;
        }
        i++;
      }else{
        sscanf(tmp,"%d %d %lf", &row, &col, &a);
        if(row == col){//diagonal elements
          nnonzero_row[row-1]++;
          i++;
        }else{//others
          nnonzero_row[row-1]++;
          nnonzero_row[col-1]++;
          i = i+2;
        }
      }
    }
  }
  nnonzero = i-1;

  row_ptr=(int *)malloc(sizeof(int)*(n+1));
  row_ptr[0] = 0;
  for(i=1;i<n+1;i++){
    row_ptr[i] = row_ptr[i-1] + nnonzero_row[i-1];
  }

  free(nnonzero_row);
  fclose(fp); //close file

  //next scan
  //file openerror
  if( (fp=fopen(argv[1],"r")) == NULL){
    printf("file open error!\n");
    exit(1);
  }

  //read file
  i=0;
  while( fgets(tmp,sizeof(tmp),fp) ){
    if(tmp[0] == '%'){
      //ignore commments
    }else{
      if(i == 0){
        sscanf(tmp,"%d %d %d", &n, &n, &j);//number of row,col,nonzero
        printf("n:%d nnonzero:%d\n", n, nnonzero);
        val = (double *)malloc(sizeof(double)*nnonzero);
        col_ind = (int *)malloc(sizeof(int)*nnonzero);
        fill = (int *)malloc(sizeof(int)*(n+1));
        for(j=0;j<nnonzero;j++){//initialize
          val[j] = 0.0;
          col_ind[j] = 0;
        }
        for(j=0;j<n+1;j++){
          fill[j] = 0;
        }
        i++;
      }else{
        sscanf(tmp,"%d %d %lf", &row, &col, &a);
        row=row-1;
        col=col-1;
        if(row != col){
          col_ind[ row_ptr[col] + fill[col] ] = row;
          val[ row_ptr[col] + fill[col] ] = a;
          fill[col]++;

          col_ind[ row_ptr[row] + fill[row] ] = col;
          val[ row_ptr[row] + fill[row] ] = a;
          fill[row]++;
        }else{
          col_ind[ row_ptr[row] + fill[row] ] = col;
          val[ row_ptr[row] + fill[row] ] = a;
          fill[row]++;
        }//endrow == col
      }//end scan row,col,a

    }//tmp[0] != 0
  }//end while

  free(fill);
  fclose(fp);

  // iccg
  int nitecg = 5000;
  double *b,*solx;
  b = (double *)malloc(sizeof(double)*n);
  solx = (double *)malloc(sizeof(double)*n);
  for(i=0;i<n;i++){
    b[i] = 1.0;
  }

  double err = 0.001;

  int *iuhead,*iucol;
  iuhead = (int *)malloc(sizeof(int)*(n+1));
  iucol = (int *)malloc(sizeof(int)*nnonzero);

  double *u;
  u = (double *)malloc(sizeof(double)*nnonzero);
  double *p,*q,*pn,*r;
  p = (double *)malloc(sizeof(double)*n);
  q = (double *)malloc(sizeof(double)*n);
  pn = (double *)malloc(sizeof(double)*n);
  r = (double *)malloc(sizeof(double)*n);
  double cgropp,cgrop;
  double alpha,alphat,beta,ar0;
  double *diag,*z;
  diag = (double *)malloc(sizeof(double)*n);
  z = (double *)malloc(sizeof(double)*n);

  int k,ku,kl,jj,ji,jjp,jp;
  double ganma,rnorm,bnorm;
  int ite;
  // convergence check
  int h,it,l,lmax,m;
  double *_solx;
  int *check;
  h = 1;
  m = 10;
  _solx = (double *)malloc(sizeof(double)*(n*m));
  check = (int *)malloc(sizeof(int)*(m));
  lmax = ceil( log(nitecg) / log(m) );

  // init
  ganma = 1.0;
  for(i=0;i<n;i++){
    diag[i] = 0.0;
    iuhead[i] = 0.0;
  }
  for(i=0;i<nnonzero;i++){
    iucol[i] = 0;
    u[i] = 0.0;
  }

  ku = 0;
  iuhead[0] = 1;

  for(i=0;i<n;i++){
    ku = 0;
    kl = 0;
    for(j=row_ptr[i];j<row_ptr[i+1];j++){
      jj = col_ind[j];
      if(jj == i){
        diag[i] = val[j] * ganma;
      }else if(jj < i){
        kl = kl + 1;
      }else if(jj>i){
        iucol[ku + row_ptr[i]] = jj;
        u[ku + iuhead[i]] = val[j];
        ku = ku + 1;
      }else{
        printf("error");
        break;
      }
      iuhead[i+1] = iuhead[i] + ku;
    }
  }


  for(i=0;i<n;i++){
    for(j=iuhead[i];j<iuhead[i+1];j++){
      jj = iucol[j];
      diag[iucol[j]] = diag[iucol[j]] - u[j] * u[j] / diag[i];

      for(jp=j+1;jp<iuhead[i+1];jp++){
        jjp = iucol[jp];
        for(ji=iuhead[jj];ji<iuhead[jj+1];ji++){
          if(iucol[ji] == jjp){
            u[ji] = u[ji] - u[j] * u[jp] / diag[i];
            //exit??
            break;
          }
        }

      }

    }
    if(abs(diag[i]) < 0.001){
      printf("diag error %d %lf\n",i,diag[i]);
      break;
    }
  }

  bnorm = 0.0;
  for(i=0;i<n;i++){
    bnorm = bnorm + abs(b[i]) * abs(b[i]);
  }

  printf("bnorm = %f , %d\n",bnorm,n);
  for(i=0;i<n;i++){
    solx[i] = 0.0;
    diag[i] = 1.0 / diag[i];
  }

  for(i=0;i<n;i++){
    ar0 = 0.0;
    for(j=row_ptr[i];j<row_ptr[i+1];j++){
      jj = col_ind[j];
      ar0 = ar0 + val[j] * solx[jj];
    }
    r[i] = b[i] - ar0;
  }

  cgrop = 0.0;

  for(ite=1;ite<nitecg;ite++){
    for(i=0;i<n;i++){
      z[i] = r[i];
    }

    for(i=0;i<n;i++){
      for(j=iuhead[i];j<iuhead[i+1];j++){
        jj = iucol[j];
        z[jj] = z[jj] - z[i] * u[j] * diag[i];
      }
    }

    z[n-1] = z[n-1] * diag[n-1];

    for(i=n-2;i>=0;--i){
      for(j=iuhead[i+1]-1;j>=iuhead[i];--j){
        jj = iucol[j];
        z[i] = z[i] - u[j] * z[jj];
      }
      z[i] = z[i] * diag[i];
    }

    for(i=0;i<n;i++){
      z[i] = r[i];
    }

    cgropp = cgrop;
    cgrop = 0.0;

    for(i=0;i<n;i++){
      cgrop = cgrop + r[i] * z[i];
    }

    if(ite==1){
      for(i=0;i<n;i++){
        pn[i] = z[i];
      }
    }else{
      beta = cgrop / cgropp;
      for(i=0;i<n;i++){
        pn[i] = z[i] + beta * p[i];
      }
    }

    for(i=0;i<n;i++){
      q[i] = 0.0;
    }

    for(i=0;i<n;i++){
      for(j=row_ptr[i];j<row_ptr[i+1];j++){
        jj = col_ind[j];
        q[i] = q[i] + val[j] * pn[jj];
      }
    }


    alphat = 0.0;
    for(i=0;i<n;i++){
      alphat = alphat + pn[i] * q[i];
    }
    alpha = cgrop / alphat;

    for(i=0;i<n;i++){
      solx[i] = solx[i] + alpha * pn[i];
      r[i] = r[i] - alpha * q[i];
    }

    for(i=0;i<n;i++){
      p[i] = pn[i];
    }

    rnorm = 0.0;

    for(i=0;i<n;i++){
      rnorm = rnorm + abs(r[i]) * abs(r[i]);
    }

    printf("ite:%d, %lf\n",ite,sqrt(rnorm / bnorm));
    if(sqrt(rnorm / bnorm) < err){
      break;
    }

    /*--- algorithm1 ---*/
    if( ite%h == 0 ){
      it = 0;
      for(l=0;l<lmax;l++){
        it = it + pow(-1,l) * floor( (ite-1) / pow(m,l) );
      }
      j = it%m;
      for(i=0;i<n;i++){
        _solx[j*n + i] = solx[i];
      }
      check[j] = ite;
      if(ite == h*m){
        h = h*2;
      }
    }
    /*--- end algorithm1 ---*/
  }//end iccg
  printf("check\n");
  for(i=0;i<m;i++){
    printf("%d\n",check[i]);
  }
  // e = x - x~
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      _solx[(i*n) + j] = solx[j] - _solx[(i*n) + j];
    }
  }
  /*--- Modified Gram-Schmidt orthogonalization ---*/
  double *enorm,*er,*eq;
  enorm = (double *)malloc(sizeof(double)*m);
  er = (double *)malloc(sizeof(double)*(m*m));
  eq = (double *)malloc(sizeof(double)*(m*n));

  // initialize enorm,er
  for(i=0;i<m;i++){
    enorm[i] = 0;
  }
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      er[i*m +j] = 0;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      enorm[i] += _solx[i*n + j] * _solx[i*n + j];
    }
    enorm[i] = sqrt(enorm[i]);
  }

  for(i=0;i<m;i++){

    er[i*m + i] = enorm[i];
    for(j=0;j<n;j++){
      eq[i*n + j] = _solx[i*n + j]/er[i*m + i];
    }

    for(j=i+1;j<m;j++){
      for(k=0;k<n;k++){
        er[i*m + j] += eq[i*n + k] * _solx[j*n +k];
      }
      for(k=0;k<n;k++){
        _solx[j*n + k] = _solx[j*n + k] - eq[i*n + k]*er[i*m + j];
      }
    }
  }
  /*--- end Modified Gram-Schmidt orthogonalization ---*/

  /*--- E^T*A*E---*/
  double *ae,*X;
  ae  = (double *)malloc(sizeof(double)*(n*m));
  X = (double *)malloc(sizeof(double)*(m*m));
  X2 = (double *)malloc(sizeof(double)*(m*m);  
//*--- init ae ---*//
  for(i=0;i<m;i++){
     for(j=0;j<n;j++){
       ae[i*n+j] = 0;
     }
   }
   for(i=0;i<m;i++){
     for(j=row_ptr[i];j<row_ptr[i+1];j++){
       jj = col_ind[j];
       ae[i*n+jj] = val[j] * eq[i*n+jj];
     }
   }
  //*--- end init ae ---*//
  // //X = eq^T * ae
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, eq, n, ae, n, 0.0, X, m);
  //E^T*ae

  /*--- end E^T*A*E---*/

  printf("# of ite. = %d, %lf",ite,sqrt(rnorm / bnorm));

  free(iuhead);
  free(iucol);
  free(u);
  free(p);
  free(q);
  free(pn);
  free(r);
  free(diag);
  free(z);
  free(_solx);
  
  free(enorm);
  free(er);
  free(eq);
  free(ae);
}
