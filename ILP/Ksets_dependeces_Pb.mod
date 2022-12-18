/*********************************************
 * OPL 12.7.0.0 Model
 * Author: Asus
 * Creation Date: 8 août 2022 at 19:03:07
 *********************************************/
 
//Paramètres et décalaration des données
int k=...; //K ensembles
range a= 1..k;
int mat_incidence[a][a]=...;
int na[a]=...; //chaque ensemble contient na objets
//int m=...; //nombre d'objet
int m=...; //le nombre d'objet maximum dans na
//range ind= 1..na[a]; 
//range ind= 1..m;
int p[a][1..m]=...;//profit de l'objet a dans l'ensemble k
int w[a][1..m]=...;//Poids de l'objet a dans l'ensemble k

int Da[a][a];//Matrice d'incidence en 0,1 la matrice égale 1 si l'ensemble i dépend directement
                             // de l'ensemble j; 0 sinon 
                             
int Ua[a][a];//Matrice d'incidence en 0,1 la matrice égale 1 si l'ensemble i dépend indirectement
                             // de l'ensemble j; 0 sinon 
                             
int W=...; //La borne des poids.
int compteur=0; 
//extraire Da du graphe en entrée  
execute{

    for(var j=1; j<=k; j++)
    {
       for(var i=1; i<=k; i++)
       {
          if(mat_incidence[i][j]==1)
          {
               Da[j][i]=1;      
          }
          else
          {
               Da[j][i]=0;      
          }
       }
    }
}

//extraire Ua du graphe en entrée  //Ua[a][a]=mat_incidence[a][a];
 
  
execute
{
    for(var j=1; j<=k;j++)
    {
        for(var i=1; i<k; i++)
        {
           if(mat_incidence[i][j]==1)
           {
               Ua[i][j]=1;         
           }
           else
           {
               Ua[i][j]=0;           
           }        
        }    
    }
         
   for(var i=1; i<=k; i++)
    {
        for(var j=1; j<=k; j++)
        {
            if((Ua[i][j]==1)&&(j!=i))
            {  compteur= compteur+1;
                 for(var b=1; b<=k; b++)
                 {
                      if((b!=j)&&(mat_incidence[j][b]==1))
                      {
                           Ua[i][b]=1;                      
                      }                 
                 }            
            }               
        }
        if(compteur ==0)
        {
            for(var j=1; j<=k;j++)
            {
                 if((Ua[i][j]==0)&&(Ua[j][i]==1))
                 {
                       Ua[i][j]=Ua[j][i];                 
                 }         
            }        
        }  
        compteur=0;  
    }
     for(var i=1; i<=k; i++)
    {
        for(var j=1; j<=k; j++)
        {
            if((Ua[i][j]==1)||(j==i))
            { writeln("i,j "+i+" "+j ); 
                for(var b=1; b<=k; b++)
                {
                 if((Ua[i][b]==0)&&(b!=i)&&(b!=j)){
                 writeln("i,j,b (e1) "+i+" "+j+" "+b );                  
                 	//writeln (mat_incidence[b][j]);
                    if((mat_incidence[b][j]==1))
                    {writeln("i,j,b "+i+" "+j+" "+b ); 
                     
                        Ua[i][b]=1; j=1;   
                       // writeln("i,j,b "+i+" "+j+" "+b );                 
                    }  
                  }                                  
                }          
            }        
        }    
    }
    
    
    
    /*for(var i=1; i<=k; i++)
    {
        for(var j=1; j<=k; j++)
        {
            if((Ua[i][j]==1)||(j==i))
            { writeln("i,j "+i+" "+j ); 
                for(var b=1; b<=k; b++)
                {
                 if((Ua[i][b]==0)&&(b!=i)&&(b!=j)){
                 writeln("i,j,b (e1) "+i+" "+j+" "+b );                  
                 	//writeln (mat_incidence[b][j]);
                    if((mat_incidence[b][j]==1))
                    {writeln("i,j,b "+i+" "+j+" "+b ); 
                     
                        Ua[i][b]=1;   
                       // writeln("i,j,b "+i+" "+j+" "+b );                 
                    }  
                  }                                  
                }          
            }        
        }    
    }*/
}


//Variables de décisions:
dvar boolean xia[a][1..m]; //la variable binaire xia pour i=1 à na et j=1 à k
dvar boolean ya[a]; // la variable ya binaire pour a allant de 1 à k
dvar int q;
//Fonction objectif

dexpr int z=sum(j in 1..k , i in 1..m) p[j][i] * xia[j][i];
maximize z;

//Contraintes:

subject to
{

q==sum (j in 1..k , i in 1..m) w[j][i]*xia[j][i] ;
sum (j in 1..k , i in 1..m) w[j][i]*xia[j][i] <=W;

//forall (j in 1..k , b in 1..k: j!=b) if(Da[j][b]==1){ya[j]-ya[b] <=0;}
forall (j in 1..k , b in 1..k: j!=b) ya[j]-ya[b] <= 1-Da[j][b];

//forall (j in 1..k , b in 1..k: j!=b) if(Ua[j][b]==0){ya[j]+ya[b] <=1;}
forall (j in 1..k , b in 1..k: j!=b) ya[j]+ya[b] <=1+Ua[j][b];

forall (j in 1..k , i in 1..na[j]) xia[j][i]-ya[j] <=0;
forall(i in 1..k)   ya[i]   <= sum(j in 1..m)xia[i][j];
forall(j in a, i in na[j]+1..m) xia[j][i]==0;
//forall(i in a) sum(j in a) Ua[i][j]*ya[j]>=1;
//xia[5][1]==1;
//ya[5]==1;  
}
