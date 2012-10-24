#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<map>
#include<sstream>
#include<pthread.h>
#include<sys/time.h>
#include<time.h>
#include<papi.h>
#define TSIZE 4

struct timeval start,end;
long seconds,useconds,mtime;


using std::copy;
using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::map;
using std::pair;
using std::stringstream;

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

long calculateTime()
{
  seconds=end.tv_sec-start.tv_sec;
  useconds=end.tv_usec-start.tv_usec;
  return ((seconds)*1000+useconds/1000.0)+0.5;
}



struct StringArray
{

    int index;
    int startindex;
    int threadid;
};

 string GeneString[100000];

map<string,string> lookup;
map<string,string>::iterator lk;

map<string,string> lookup2;
map<string,string>::iterator lk2;

map<string,string> lookup3;
map<string,string>::iterator lk3;

map<string,string> lookup4;
map<string,string>::iterator lk4;




map<char,int> gene;
map<char,int>::iterator ge;



int EditMatrix[1000][1000];
char tempDNA[40];
string DNAStringArray;
char tempNumbers[40];
char InputString_1[100]={'\0'};
char InputString_2[100]={'\0'};
int InputCount1=0,InputCount2=0;


vector<string> NumbersStringArray(500000);

const unsigned long num_elements_Numbers = 3;
const unsigned long num_elements_DNA = 4;
int comb_len;
int Rsolution[10],Csolution[10];
int ncount=0;


void Permute_DNA(char *,char *);
void Permute_Numbers(char *,char *);
void *DivideTask(void *);
void EditDistance1(string );
void Dictionary1(string,string);
void EditDistance2(string );
void Dictionary2(string,string);
void EditDistance3(string );
void Dictionary3(string,string);
void EditDistance4(string );
void Dictionary4(string,string);
void CombineDictionary();
void *Distribute(void *);
void Fillup(int ,int,string);
void test();

void CombinationsRecursiveNumbers(const vector<char> &elems, unsigned long req_len,
			    vector<unsigned long> &pos, unsigned long depth,
			    unsigned long margin)
{

	if (depth >= req_len) {
		for (unsigned long ii = 0; ii < pos.size(); ++ii)
			{


			tempNumbers[ii]=elems[pos[ii]];

			}
       Permute_Numbers(tempNumbers,tempNumbers);

		return;
	}


	for (unsigned long ii = margin; ii < elems.size(); ++ii) {
		pos[depth] =ii;
		CombinationsRecursiveNumbers(elems, req_len, pos, depth + 1, ii);
	}
	return;
}


void Permute_Numbers(char *string_start, char *p) {

  if (*(p+1) == 0)
  {
     NumbersStringArray[ncount]=string_start;
     ncount++;
   }

  else {
   char *swap;

    for(swap = p; *swap; ++swap) {
      char *same;

      for(same = p; *same != *swap; ++same) {
      }
      if (same == swap) {
	char tmp = *swap;
	*swap = *p;
	*p = tmp;
	Permute_Numbers(string_start, p+1);
	*p = *swap;
	*swap = tmp;
      }
    }
  }
}


void combinations_r_numbers(const vector<char> &elems, unsigned long req_len)
{

	vector<unsigned long> positions(req_len, 0);
	CombinationsRecursiveNumbers(elems, req_len, positions, 0, 0);
}

int dcount=0;

void CombinationsRecursiveDNA(const vector<char> &elems, unsigned long req_len,
			    vector<unsigned long> &pos, unsigned long depth,
			    unsigned long margin)
{

	if (depth >= req_len) {
		for (unsigned long ii = 0; ii < pos.size(); ++ii)
			{


			tempDNA[ii]=elems[pos[ii]];

			}
       Permute_DNA(tempDNA,tempDNA);

		return;
	}


	for (unsigned long ii = margin; ii < elems.size(); ++ii) {
		pos[depth] =ii;
		CombinationsRecursiveDNA(elems, req_len, pos, depth + 1, ii);
	}
	return;
}

int maincount=0;
void Permute_DNA(char *string_start, char *p) {

int numbercount=ncount;


struct StringArray *SAA;

SAA=(struct StringArray *)malloc(sizeof(struct StringArray)*9);



int k=0;



  if (*(p+1) == 0) {
     int i=0;

     maincount=0;



     while(i<numbercount)
     {
         DNAStringArray="";
         DNAStringArray=string_start;
         DNAStringArray=DNAStringArray+NumbersStringArray[i];
         GeneString[i]=DNAStringArray;

         i++;
         maincount++;


     }

     int additional=0;
	int blocksize=0;
	int arrayoffset=0;
	int jj=0;

    pthread_t thread[TSIZE];


	blocksize=maincount/TSIZE;



    additional=maincount%TSIZE;




	for(k=0; k<TSIZE; k++)
	{
	     if(k<additional)
	     {
	         SAA[k].index=blocksize+1;


	     }
	     else
	     {
	         SAA[k].index=blocksize;

	     }

	     SAA[k].startindex=arrayoffset;
	     SAA[k].threadid=k;
	     pthread_create(&thread[k], NULL ,Distribute,(void *)(SAA+k));
	     arrayoffset=arrayoffset+SAA[k].index;

	}

     for(k=0; k< TSIZE; k++)
     {
        pthread_join(thread[k], NULL);
     }
}



  else {
   char *swap;

    for(swap = p; *swap; ++swap) {
      char *same;

      for(same = p; *same != *swap; ++same) {
      }
      if (same == swap) {
	char tmp = *swap;
	*swap = *p;
	*p = tmp;
	Permute_DNA(string_start, p+1);
	*p = *swap;
	*swap = tmp;
      }
    }
  }



}


void test()
{
    int t1=0;

    int sum=0;

    for ( lk=lookup.begin() ; lk!=lookup.end(); lk++ )
      {
         // cout<<endl;
         // cout<<(*lk).first<<" ||"<<(*lk).second;
          t1++;
       }

       sum=t1;
//cout<<" "<<endl;
cout<<"\nThe Size of The Dictionary for The T-Block- "<<comb_len<<" :"<<sum<<endl;
}

int INDEX=0;

void *Distribute(void *arg)
{
   int loop1=0;
   int loop2=0;
   int loop3=0;
   int loop4=0;
   int loop5=0;
   int loop6=0;
   int loop7=0;
   int loop8=0;


   struct StringArray *SAPTR = (struct StringArray *)arg;


    if(SAPTR->threadid==0)
    {
        for(loop1=SAPTR->startindex; loop1<SAPTR->startindex+SAPTR->index; loop1++)
        {
            pthread_mutex_lock(&lock);

            EditDistance1(GeneString[loop1]);
            pthread_mutex_unlock(&lock);
        }

    }
    if(SAPTR->threadid==1)
    {
        for(loop2=SAPTR->startindex; loop2<SAPTR->startindex+SAPTR->index; loop2++)
        {
            pthread_mutex_lock(&lock);
            EditDistance2(GeneString[loop2]);
            pthread_mutex_unlock(&lock);
        }

    }
    if(SAPTR->threadid==2)
    {
        for(loop3=SAPTR->startindex; loop3<SAPTR->startindex+SAPTR->index; loop3++)
        {
            pthread_mutex_lock(&lock);
            EditDistance3(GeneString[loop3]);
            pthread_mutex_unlock(&lock);
        }

    }
    if(SAPTR->threadid==3)
    {
        for(loop4=SAPTR->startindex; loop4<SAPTR->startindex+SAPTR->index; loop4++)
        {
            pthread_mutex_lock(&lock);
            EditDistance4(GeneString[loop4]);
            pthread_mutex_unlock(&lock);
        }


    }
}

void EditDistance1(string String)
{
int loop=0;

    int length=comb_len;
    char* input = ( char*)malloc( sizeof( char ) *((String).length() +1) );
    strcpy( input,(String).c_str() );

    int EditArray[200][200];

    int i=0,j=0;

for(i=0;i<length;i++)
{
    for(j=0;j<length;j++)
    {
        EditArray[i][j]=0;
    }
}

char Row[20],Column[20];
for(i=0;i<length-1;i++)
 {
     Row[i]=input[i];
 }

 j=0;

 i=0;

 for(i=length-1;i<2*(length-1);i++)
 {
     Column[j]=input[i];
     j++;
 }

 j=1;
 i=0;
int check=0;
 for(i=2*(length-1);i<3*(length-1);i++)
 {
      char a =input[i];
      check=0;
      check=atoi(&a);
      if(check==3)
      {
        check=-1;
        EditArray[j][0]=EditArray[j-1][0]+check;
      }
      else
      {
        EditArray[j][0]=EditArray[j-1][0]+check;
      }

     j++;
 }
 j=1;
 i=0;

 for(i=3*(length-1);i<4*(length-1);i++)
 {
      char b =input[i];
      check=0;
      check=atoi(&b);
      if(check==3)
      {
        check=-1;
        EditArray[0][j]=EditArray[0][j-1]+check;
      }
      else
      {
        EditArray[0][j]=EditArray[0][j-1]+check;
      }

     j++;

 }

 for (i=1;i<length;i++)
   {
       for(j=1;j<length;j++)
        {
            if(Row[i-1]==Column[j-1])
            {
                EditArray[i][j]=EditArray[i-1][j-1];
            }
            else{
                  EditArray[i][j]=min(EditArray[i][j-1]+1,min(EditArray[i-1][j-1]+1,EditArray[i-1][j]+1));
                }
        }
   }

   j=0;

   for(i=0;i<length-1;i++)
   {

       Csolution[j]=(EditArray[length-1][i+1])-(EditArray[length-1][i]);
       Rsolution[j]=(EditArray[i+1][length-1])-(EditArray[i][length-1]);
       j++;
   }


   string tempstring="";



for(i=0;i<length-1;i++)
{
    stringstream ss;
    ss<<Rsolution[i];
    tempstring=tempstring+ss.str();
}
for(i=0;i<length-1;i++)
{
    stringstream ss1;
    ss1<<Csolution[i];
    tempstring=tempstring+ss1.str();
}

Dictionary1(String,tempstring);





}


void Dictionary1(string s,string x)
 {

  lk=lookup.begin();
  lookup.insert(lk,pair<string,string>(s,x));

 }


void EditDistance2(string String)
{
int loop=0;

    int length=comb_len;
    char* input = ( char*)malloc( sizeof( char ) *((String).length() +1) );
    strcpy( input,(String).c_str() );

    int EditArray[200][200];

    int i=0,j=0;

for(i=0;i<length;i++)
{
    for(j=0;j<length;j++)
    {
        EditArray[i][j]=0;
    }
}

char Row[20],Column[20];
for(i=0;i<length-1;i++)
 {
     Row[i]=input[i];
 }

 j=0;

 i=0;

 for(i=length-1;i<2*(length-1);i++)
 {
     Column[j]=input[i];
     j++;
 }

 j=1;
 i=0;
int check=0;
 for(i=2*(length-1);i<3*(length-1);i++)
 {
      char a =input[i];
      check=0;
      check=atoi(&a);
      if(check==3)
      {
        check=-1;
        EditArray[j][0]=EditArray[j-1][0]+check;
      }
      else
      {
        EditArray[j][0]=EditArray[j-1][0]+check;
      }

     j++;
 }
 j=1;
 i=0;

 for(i=3*(length-1);i<4*(length-1);i++)
 {
      char b =input[i];
      check=0;
      check=atoi(&b);
      if(check==3)
      {
        check=-1;
        EditArray[0][j]=EditArray[0][j-1]+check;
      }
      else
      {
        EditArray[0][j]=EditArray[0][j-1]+check;
      }

     j++;

 }

 for (i=1;i<length;i++)
   {
       for(j=1;j<length;j++)
        {
            if(Row[i-1]==Column[j-1])
            {
                EditArray[i][j]=EditArray[i-1][j-1];
            }
            else{
                  EditArray[i][j]=min(EditArray[i][j-1]+1,min(EditArray[i-1][j-1]+1,EditArray[i-1][j]+1));
                }
        }
   }

   j=0;

   for(i=0;i<length-1;i++)
   {

       Csolution[j]=(EditArray[length-1][i+1])-(EditArray[length-1][i]);
       Rsolution[j]=(EditArray[i+1][length-1])-(EditArray[i][length-1]);
       j++;
   }


   string tempstring="";



for(i=0;i<length-1;i++)
{
    stringstream ss;
    ss<<Rsolution[i];
    tempstring=tempstring+ss.str();
}
for(i=0;i<length-1;i++)
{
    stringstream ss1;
    ss1<<Csolution[i];
    tempstring=tempstring+ss1.str();
}

Dictionary2(String,tempstring);





}

void Dictionary2(string s,string x)
 {

  lk2=lookup2.begin();
  lookup2.insert(lk2,pair<string,string>(s,x));
 }


void EditDistance3(string String)
{
int loop=0;

    int length=comb_len;
    char* input = ( char*)malloc( sizeof( char ) *((String).length() +1) );
    strcpy( input,(String).c_str() );

    int EditArray[200][200];

    int i=0,j=0;

for(i=0;i<length;i++)
{
    for(j=0;j<length;j++)
    {
        EditArray[i][j]=0;
    }
}

char Row[20],Column[20];
for(i=0;i<length-1;i++)
 {
     Row[i]=input[i];
 }

 j=0;

 i=0;

 for(i=length-1;i<2*(length-1);i++)
 {
     Column[j]=input[i];
     j++;
 }

 j=1;
 i=0;
int check=0;
 for(i=2*(length-1);i<3*(length-1);i++)
 {
      char a =input[i];
      check=0;
      check=atoi(&a);
      if(check==3)
      {
        check=-1;
        EditArray[j][0]=EditArray[j-1][0]+check;
      }
      else
      {
        EditArray[j][0]=EditArray[j-1][0]+check;
      }

     j++;
 }
 j=1;
 i=0;

 for(i=3*(length-1);i<4*(length-1);i++)
 {
      char b =input[i];
      check=0;
      check=atoi(&b);
      if(check==3)
      {
        check=-1;
        EditArray[0][j]=EditArray[0][j-1]+check;
      }
      else
      {
        EditArray[0][j]=EditArray[0][j-1]+check;
      }

     j++;

 }

 for (i=1;i<length;i++)
   {
       for(j=1;j<length;j++)
        {
            if(Row[i-1]==Column[j-1])
            {
                EditArray[i][j]=EditArray[i-1][j-1];
            }
            else{
                  EditArray[i][j]=min(EditArray[i][j-1]+1,min(EditArray[i-1][j-1]+1,EditArray[i-1][j]+1));
                }
        }
   }

   j=0;

   for(i=0;i<length-1;i++)
   {

       Csolution[j]=(EditArray[length-1][i+1])-(EditArray[length-1][i]);
       Rsolution[j]=(EditArray[i+1][length-1])-(EditArray[i][length-1]);
       j++;
   }


   string tempstring="";



for(i=0;i<length-1;i++)
{
    stringstream ss;
    ss<<Rsolution[i];
    tempstring=tempstring+ss.str();
}
for(i=0;i<length-1;i++)
{
    stringstream ss1;
    ss1<<Csolution[i];
    tempstring=tempstring+ss1.str();
}

Dictionary3(String,tempstring);





}


void Dictionary3(string s,string x)
 {

  lk3=lookup3.begin();
  lookup3.insert(lk3,pair<string,string>(s,x));
 }


 void EditDistance4(string String)
{
int loop=0;

    int length=comb_len;
    char* input = ( char*)malloc( sizeof( char ) *((String).length() +1) );
    strcpy( input,(String).c_str() );

    int EditArray[200][200];

    int i=0,j=0;

for(i=0;i<length;i++)
{
    for(j=0;j<length;j++)
    {
        EditArray[i][j]=0;
    }
}

char Row[20],Column[20];
for(i=0;i<length-1;i++)
 {
     Row[i]=input[i];
 }

 j=0;

 i=0;

 for(i=length-1;i<2*(length-1);i++)
 {
     Column[j]=input[i];
     j++;
 }

 j=1;
 i=0;
int check=0;
 for(i=2*(length-1);i<3*(length-1);i++)
 {
      char a =input[i];
      check=0;
      check=atoi(&a);
      if(check==3)
      {
        check=-1;
        EditArray[j][0]=EditArray[j-1][0]+check;
      }
      else
      {
        EditArray[j][0]=EditArray[j-1][0]+check;
      }

     j++;
 }
 j=1;
 i=0;

 for(i=3*(length-1);i<4*(length-1);i++)
 {
      char b =input[i];
      check=0;
      check=atoi(&b);
      if(check==3)
      {
        check=-1;
        EditArray[0][j]=EditArray[0][j-1]+check;
      }
      else
      {
        EditArray[0][j]=EditArray[0][j-1]+check;
      }

     j++;

 }

 for (i=1;i<length;i++)
   {
       for(j=1;j<length;j++)
        {
            if(Row[i-1]==Column[j-1])
            {
                EditArray[i][j]=EditArray[i-1][j-1];
            }
            else{
                  EditArray[i][j]=min(EditArray[i][j-1]+1,min(EditArray[i-1][j-1]+1,EditArray[i-1][j]+1));
                }
        }
   }

   j=0;

   for(i=0;i<length-1;i++)
   {

       Csolution[j]=(EditArray[length-1][i+1])-(EditArray[length-1][i]);
       Rsolution[j]=(EditArray[i+1][length-1])-(EditArray[i][length-1]);
       j++;
   }


   string tempstring="";



for(i=0;i<length-1;i++)
{
    stringstream ss;
    ss<<Rsolution[i];
    tempstring=tempstring+ss.str();
}
for(i=0;i<length-1;i++)
{
    stringstream ss1;
    ss1<<Csolution[i];
    tempstring=tempstring+ss1.str();
}

Dictionary4(String,tempstring);




}


void Dictionary4(string s,string x)
 {

  lk4=lookup4.begin();
  lookup4.insert(lk4,pair<string,string>(s,x));
 }

void CombineDictionary()
{
    for ( lk2=lookup2.begin() ; lk2 !=lookup2.end(); lk2++ )
      {
           string temp1=(*lk2).first;
           if(temp1==(*lk2).first)
           {
           string temp2=(*lk2).second;
           lk=lookup.begin();
           lookup.insert(lk,pair<string,string>(temp1,temp2));
           }
      }
    for ( lk3=lookup3.begin() ; lk3 !=lookup3.end(); lk3++ )
      {
           string temp3=(*lk3).first;
           if(temp3==(*lk3).first)
           {

           string temp4=(*lk3).second;
           lk=lookup.begin();
           lookup.insert(lk,pair<string,string>(temp3,temp4));
           }
      }
    for ( lk4=lookup4.begin() ; lk4 !=lookup4.end(); lk4++ )
      {
           string temp5=(*lk4).first;
           if(temp5==(*lk4).first)
           {

           string temp6=(*lk4).second;
           lk=lookup.begin();
           lookup.insert(lk,pair<string,string>(temp5,temp6));
           }
      }

}



void GeneMapping(char s,int x)
{
  ge=gene.begin();
gene.insert(ge,pair<char,int>(s,x));
}

int PrevColumn[100]={0};
int PrevCount=0;

void FourRussian()
{
	int i=0, j=0, k=0;


        cout<<""<<endl;
        for(i = 0; (InputString_2[i]=getchar())!='\n';)
            {
            }

         cout<<"Enter the Input String 1:"<<endl;
        for(i = 0; (InputString_1[i]=getchar())!='\n';)
            {
                 for ( ge=gene.begin() ; ge !=gene.end(); ge++ )
	             {
	                  if((*ge).first==InputString_1[i])
	                  {
	                      EditMatrix[i+2][0]=(*ge).second; //ROW INPUT
	                  }
                 }
                 i++;

            }

        InputString_1[i] = '\0';
        InputCount1=i;

        cout<<"Enter the Input String 2:"<<endl;
        for(i = 0; (InputString_2[i]=getchar())!='\n';)
            {
                 for ( ge=gene.begin() ; ge !=gene.end(); ge++ )
	             {
	                  if((*ge).first==InputString_2[i])
	                  {
	                      EditMatrix[0][i+2]=(*ge).second;   //COLUMN INPUT
	                  }
                 }
                 i++;
             }
             InputString_2[i] = '\0';
             InputCount2=i;


        EditMatrix[0][1]=-9;
        EditMatrix[1][0]=-9;

     for(i=1;i<InputCount1+2;i++)
     {
         EditMatrix[i][1]=1;
         EditMatrix[1][i]=1;

     }

     EditMatrix[0][0]=0;
     EditMatrix[1][1]=0;




   int T=comb_len;
   int VSlider=1;
   int HSlider=1;


int OverwriteHSlider=0;
int PrevVSliderCount=0;
int StartingPoint=0;

  while(VSlider<=InputCount1+1)
  {

   if(HSlider>=T && PrevVSliderCount>=T)
   {
       OverwriteHSlider=T;

       int tempcount=0;

       for(StartingPoint=OverwriteHSlider;StartingPoint<=InputCount1+1;StartingPoint+=(OverwriteHSlider-1))

       {
          EditMatrix[PrevVSliderCount][StartingPoint]=PrevColumn[StartingPoint];

        }

       PrevVSliderCount=0;
       PrevCount=0;
   }



    HSlider=1;

   while(HSlider<=(InputCount2-1))
   {
     string match="";

     for(j=VSlider;j<VSlider+(T-1);j++)
      {
        for ( ge=gene.begin() ; ge !=gene.end(); ge++ )
         {
          if(EditMatrix[j+1][0]==(*ge).second)
          {

           match=match+(*ge).first;

          }
       }

     }

      for(j=HSlider;j<HSlider+(T-1);j++)
      {
          for ( ge=gene.begin() ; ge !=gene.end(); ge++  )
         {
          if(EditMatrix[0][j+1]==(*ge).second)
          {

           match=match+(*ge).first;

          }
        }
      }

      for(j=VSlider;j<VSlider+(T-1);j++)
      {
             if(EditMatrix[j+1][HSlider]==-1)
             {
                 stringstream rowno3;
                 rowno3<<3;
                 match=match+rowno3.str();
             }
             else{

            stringstream rowno;
            rowno<<EditMatrix[j+1][HSlider];
            match=match+rowno.str();

             }

      }
      for(j=HSlider;j<HSlider+(T-1);j++)
      {
          if(EditMatrix[VSlider][j+1]==-1)
          {
             stringstream colno3;
             colno3<<3;
             match=match+colno3.str();
          }
          else{
          stringstream colno;
             colno<<EditMatrix[VSlider][j+1];
             match=match+colno.str();
          }
      }

                int i=0;
                int PasteArray[100]={0};
                int cnt=0,cnt1=0;
                int flag=0;
                string temp1,temp2,temp3,temp4;

  for ( lk=lookup.begin() ; lk !=lookup.end(); lk++ )
      {
          if(match==(*lk).first)
          {
             string   temp=(*lk).second;
                const char *a=temp.c_str();


                Fillup(VSlider,HSlider,temp);


          }
       }

HSlider=HSlider+(T-1);

   }

  VSlider=VSlider+(T-1);
  PrevVSliderCount=VSlider;

  }


  for(i=0;i<InputCount1+2;i++)
{
    cout<<"  "<<endl;
    for(j=0;j<InputCount2+2;j++)
    cout<<" | "<<EditMatrix[i][j];
}

int columnsum=0;
int rowsum=0;

for(i=0;i<InputCount1+2;i++)
{
  rowsum+=EditMatrix[i+2][1];
}

for(i=0;i<InputCount2+2;i++)
{
    columnsum+=EditMatrix[InputCount1+1][i+2];
}
cout<<" "<<endl;

cout<<"\nThe Edit Distance of the NXN Matrix is  :"<<columnsum+rowsum<<endl;
}

void Fillup(int vslider,int hslider,string temp)
{

    const char *a=temp.c_str();

    int cnt=0;
    int fillcount=0;
    int PasteArray[100]={2};
    int i=0,j=0;

                while(a[cnt]!='\0')
                {

                    if(a[cnt]=='-')
                    {
                       cnt++;

                       PasteArray[fillcount]=-1;

                      cnt++;

                       fillcount++;
                    }
                    else
                    {
                     char atemp=a[cnt];

                     PasteArray[fillcount]=atoi(&atemp);

                     cnt++;
                     fillcount++;
                    }
                 }
int T=comb_len;
int k=0;
int storerow=0;
int storecol=0;

   for(i=vslider;i<vslider+(T-1);i++)
     {


         EditMatrix[i+1][hslider+(T-1)]=PasteArray[k];
         storerow=i+1;
         storecol=hslider+(T-1);
         k++;

    }

   j=k;

    for(i=hslider;i<hslider+(T-1);i++)
   {
        EditMatrix[vslider+(T-1)][i+1]=PasteArray[j];
        PrevColumn[i+1]=EditMatrix[vslider+(T-1)][i+1];
        j++;
       }
   EditMatrix[storerow][storecol]=PasteArray[k-1];


}


void combinations_r(const vector<char> &elems, unsigned long req_len)
{
	vector<unsigned long> positions(req_len, 0);
	CombinationsRecursiveDNA(elems, req_len, positions, 0, 0);
}



int main()
{
    int s,temp=4;

    cout<<"Enter the T block size :"<<endl;
    cin>>comb_len;

    gettimeofday(&start,NULL);
	float rtime, ptime, ipc;

	int rval;
	long_long ins; 

	rval =PAPI_ipc(&rtime, &ptime, &ins,  &ipc );
	if(rval != PAPI_OK)
	{
		perror("PAPI_ipc error\n");
	}




	vector<char> elements_Numbers(num_elements_Numbers);
	char elements_str_Numbers[num_elements_Numbers + 1] = "103";
	copy(elements_str_Numbers, elements_str_Numbers+ num_elements_Numbers, elements_Numbers.begin());
	combinations_r_numbers(elements_Numbers, (comb_len*2)-2);

	vector<char> elements_DNA(num_elements_DNA);
	char elements_str_DNA[num_elements_DNA + 1] = "ATGC";

     for(s=0;s<4;s++)
     {
       GeneMapping(elements_str_DNA[s],temp);
       temp++;
     }



       copy(elements_str_DNA, elements_str_DNA + num_elements_DNA, elements_DNA.begin());
       combinations_r(elements_DNA, (comb_len*2)-2);
       gettimeofday(&end,NULL);
       long double time=calculateTime();
       cout<<"\n\n TIME:"<<time<<"milliseconds"<<endl;
	rval =PAPI_ipc(  &rtime, &ptime,  &ins, &ipc );
	if(rval != PAPI_OK)
	{
		perror("PAPI_ipc error\n");
	}
	printf("\nNumber of Instructions =  %lld \n", ins);
	printf("Realtime %f \n" , rtime);
	printf("Total Process Time %f \n" , ptime);
	printf("Instructions Per Cycle %f \n" , ipc);

/*
         for ( ge=gene.begin() ; ge !=gene.end(); ge++ )
         {
          cout<<endl;
          cout<<(*ge).first<<" ||"<<(*ge).second;
         }

    cout<<""<<endl;
*/
    CombineDictionary();
    test();
    FourRussian();
    cout<<" "<<endl;
    return 0;
}



