#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<iomanip>
#include<vector>
#include <typeinfo>
#include <algorithm>
using namespace std;

const int sims=1;
const int Tmax=80000;
const double equil_tim=50000; //time for an equilibration must be even, Tmax>equil_tim

const int inte_stat=1; //1 if interaction ,0 if no interaction
const int ran_wal=1; //1 if random walk ,0 if no random walk

double kbt=1; //temperature*kb
const double gamm=2; //friction coef.
double bb=5.0; //Potential strength
double diff=kbt/gamm; //diffusion coef dx*dx/dt*0.5
const double dt=1.0; // timestep
double dx=sqrt(dt*2*diff);


const int maxtemp=5.0;//Cycling temp from 1 to 5 and back
const int temp_dif=maxtemp-gamm*diff;
const double equil_tim_ha=equil_tim/2.0;
const double temp_inc=temp_dif/equil_tim_ha;


int start=0; //start=0 means equilibration and start=1 sim

int kn_nn_ori;//Number of Knot
int kn_nn; // Number of knots
const double pp=0.0001; // probability of knotting from topo
const int knots_cons=8;

double topo_le=1.0; //half of topo length
double kn_ll; //Lenght of Knot

// Interaction parameters
const double aa=0.0018;
#define pi 3.14159265
const double lb=26;
const double rb=81; //80.955
const double lzer=100.0;
const double cutoff=lzer+50.0;
double dis_kns;
double pos_kns_new1;
double pos_kns_new2;
double ran1;
double ran2;
double dis_kns_new;
double free_dis;
double add_ran;

double stds[Tmax]; //standard deviation

int ll; // Length of Polymer
double ll_ha;
int tps_size; // Number of topos

int types[knots_cons]; //Types of knots determined by their crossings

double tps;
double tt; // Time
double dummy;
int boo[Tmax]; //Booleans
int boolle[Tmax][sims];
double prob[Tmax]; // Final probabilities of beinng knotted with time
double time_arr[Tmax];
int knotted[Tmax];
int step;
double dis;
double ave;
double topo_prob[knots_cons][knots_cons];
int hit_nn;

double modmod(double num,double ll);
double theta(double x);
double ff(double x);
const double m=(ff(rb)-ff(lb))/(rb-lb);
double fr(double x);


double stand(double x[],int cut);

double probab(double (&topo_pr)[8][8],double ran,int typ);


int main(void){


  srand(time(NULL));
  for(int k=0;k<Tmax;k++) boo[k]=0; // Initialising the booleans
  for(int k=0;k<Tmax;k++) {for(int j=0;j<sims;j++) boolle[k][j]=0;}
  for(int k=0;k<Tmax;k++) stds[k]=0;

  //Open files and imopr details.txt and knot_det.txt
  ifstream stream;
  string x;
  stream.open("./details.txt");

  for(int i=0;i<2;i++) {getline(stream, x);}
  for(int i=0;i<3;i++){
    if(i==0) {
      stream>>x;
      ll=atoi(x.c_str());
    }
    if(i==1) {
      stream>>x;
      tps_size=atoi(x.c_str());
    }
    if(i==2) {
      stream>>x;
      kn_ll=atoi(x.c_str());
    }
    if(i==2) {
      stream>>x;
      kn_nn_ori=atoi(x.c_str());
    }
  }

  vector <int> type_sim_ori (kn_nn_ori);

  for(int i=0;i<kn_nn_ori;i++){
    stream>>x;
    type_sim_ori[i]=atoi(x.c_str());
  }
  stream.close();

  stream.open("./knot_det_direct_simp.txt");
  for(int i=0;i<2;i++) {getline(stream, x);}

  for(int i=0;i<knots_cons;i++){
    for(int j=0;j<knots_cons+1;j++){
      stream>>x;
      if(j==0) types[i]=atoi(x.c_str());
      if (j!=0) topo_prob[i][j-1]=atof(x.c_str());
  }
  }

  stream.close();
  ///Ending the importing of files

  // NORMALISING Probabilities/////

  double norm;
  for(int i=0;i<knots_cons;i++){
    norm=0;
    for(int j=0;j<knots_cons;j++) norm+=topo_prob[j][i];
    for(int j=0;j<knots_cons;j++) topo_prob[j][i]=topo_prob[j][i]/norm;
  }


  ////////////////////
  /////MAIN LOOP//////
  ////////////////////
  ll_ha=ll/2.0;
  double y; //Just a dummy variable to be used with probab()
  int removals;


  //Variables for free energy testing
  // double distances[Tmax];
  // double distances_sims[sims][Tmax];
  // int kkk[sims];
  // for(int i=0;i<sims;i++) kkk[i]=0;

  double pos1[Tmax];
  double pos2[Tmax];
  double dis12[Tmax];
  double dx1[Tmax];
  double dx2[Tmax];
  double cor[Tmax];
  double cor_coef=0;
  double mean1=0;
  double mean2=0;
  int aaa=0;
  for(int k=0;k<Tmax;k++) cor[k]=0;

  for(int j=0;j<sims;j++){
    int test=0;
    kn_nn=kn_nn_ori;
    vector <int> type_sim(kn_nn_ori);
    for(int i=0;i<kn_nn;i++) type_sim[i]=type_sim_ori[i];
    vector<double> cns (kn_nn);
    vector<double> displacements (kn_nn);
    vector<double> v1(2);
    vector<int> hit;
    vector< vector<double> > cros(kn_nn,v1);
    vector<int> changes;
    vector<int> remv;



    hit_nn=0;
    step=0;
    tt=0.0;
    removals=0;
    start=0;

    //Generating positions of centres of knots
    cns[0]=rand()*1.0/RAND_MAX*ll;
    for(int i=1;i<kn_nn;i++){
      cns[i]=modmod(cns[0]+ll*(1.0*i)/kn_nn,ll);
    }
    //Generating appropriate crossings at each end of the knots
    for(int i=0;i<kn_nn;i++){
      cros[i][0]=modmod(cns[i]-kn_ll/2.0,ll);
      cros[i][1]=modmod(cns[i]+kn_ll/2.0,ll);
    }
    //Generating positions of topos
    // tps=rand()*1.0/RAND_MAX*ll;
    if(cros[kn_nn-1][1]<cros[0][0]) dummy=abs(cros[kn_nn-1][1]-cros[0][0]);
    else dummy=ll-abs(cros[kn_nn-1][1]-cros[0][0]);
    tps=modmod((cros[kn_nn-1][1]+rand()*1.0/RAND_MAX*dummy),ll);


    cout<<j<<endl;


    /////////////////////////
    ///Starting simulation///
    /////////////////////////

    for(step;step<Tmax;step++){

      pos1[step]=cns[0];
      pos2[step]=cns[1];
      dis12[step]=abs(cns[1]-cns[0]);
      if(dis12[step]>ll_ha) dis12[step]=ll-dis12[step];
      if(start==1  && test==0) {
        if( dis12[step]<105){
        // cout<<  dis12[step]<<" Success"<<endl;
        aaa+=1;
        test=1;}
      }
      diff=kbt/gamm;
      dx=sqrt(dt*2*diff);;

      // if (tt>(Tmax-100) && kn_nn>1) cout <<cns[0] << " "<< cns[1]<<endl;
      if(start==0 && tt>=equil_tim) {start+=1; step=0; tt=0;}

      // if(start==1 && kn_nn>1) {
      //
      //   distances[step]=abs(cns[0]-cns[1]);
      //   if(distances[step]>ll_ha) distances[step]=ll-distances[step];
      //   distances_sims[j][step]=distances[step];
        // kkk[j]+=1;

      // }


      if(start==0){    //Used to equilibrate the system first

          if(tt<equil_tim_ha) kbt+=temp_inc;
          if(tt>=equil_tim_ha) kbt-=temp_inc;
        }


      if (j==0) time_arr[step]=tt;
      if(kn_nn>=1) knotted[step]=1;
      if(kn_nn==0) knotted[step]=0;
      tt+=dt;

      ////////////////////////////
      //Interaction-MH algorithm//
      ////////////////////////////

      if (kn_nn>=1){

        for(int i=0;i<kn_nn;i++) displacements[i]=0;

        if(inte_stat==1 && kn_nn>1){
        for(int u=0;u<(kn_nn-1);u++){
          for(int r=u+1;r<kn_nn;r++){

            dis_kns=abs(cns[u]-cns[r]);
            if (dis_kns>ll_ha) {dis_kns=ll-dis_kns;}
            if (dis_kns>cutoff) break;

            free_dis=fr(dis_kns);

            do {

              if(rand()*1.0/RAND_MAX<0.5) {
                ran1=dx;
              }
              else{
                ran1=-dx;
              }
              if(rand()*1.0/RAND_MAX<0.5) {
                ran2=dx;
              }
              else{
                ran2=-dx;
              }
              pos_kns_new1=modmod(cns[u]+ran1,ll);
              pos_kns_new2=modmod(cns[r]+ran2,ll);
              dis_kns_new=abs(pos_kns_new1-pos_kns_new2);
              if (dis_kns_new>ll_ha) {dis_kns_new=ll-dis_kns_new;}

            } while (exp((free_dis-fr(dis_kns_new))/kbt)<rand()*1.0/RAND_MAX);

              displacements[u]+=ran1;
              displacements[r]+=ran2;

            }
        }
        }
        ////////////////////////
        ///END OF INTERACTION///
        ////////////////////////

        //Moving knots randomly
        for(int i=0;i<kn_nn;i++){
          if(ran_wal==1){
          if (rand()*1.0/RAND_MAX<0.5){
            displacements[i]+=dx;
          }
          else{
            displacements[i]-=dx;
          }
        }
        }

        for(int i=0;i<kn_nn;i++) {
           cns[i]=modmod(cns[i]+displacements[i],ll);
           for(int l=0;l<2;l++) cros[i][l]=modmod(cros[i][l]+displacements[i],ll);
        }


      //Check if knots were hit by TOPO II,
      // Need to equilibrate first when start==1
      if(start==2){
      for(int i=0;i<kn_nn;i++) {
        for(int k=0;k<2;k++){
          dis=abs(tps-cros[i][k]);
          if (dis>ll_ha) dis=ll-dis;
          if(dis<topo_le) {

            hit.push_back(i);
            hit_nn+=1;

          }
        }
      }
      }
      //////////////////////
      ///END OF MOVEMENTS///
      //////////////////////

      // If any knots were hit decide whether they should change or no
      if(hit_nn!=0) {

        for(int i=0;i<hit_nn;i++){

          y=probab(topo_prob,rand()*1.0/RAND_MAX,type_sim[hit[i]]);

          if(y!=0) {

            type_sim[hit[i]]=y;

          }
          if(y==0){

            removals+=1;
            remv.push_back(hit[i]);

          }
        }
        //Remove knots if that was the case in the previous step
        if(removals!=0){

          sort(remv.begin(),remv.end());
          for(int i=0;i<removals;i++){

            kn_nn-=1;
            type_sim.erase(type_sim.begin()+remv[removals-i-1]);
            cns.erase(cns.begin()+remv[removals-i-1]);

            cros.erase(cros.begin()+remv[removals-i-1]);
            displacements.erase(displacements.begin()+remv[removals-i-1]);
          }
          removals=0;
          remv.clear();
        }

        hit.clear();
        hit_nn=0;
      }
    }

    //Generate a knot from topo when there are no knots
      // if (kn_nn==0){
      //   if (rand()*1.0/RAND_MAX<pp){
      //
      //   cns.resize(1);
      //   cns[0]=(modmod((tps+((rand()*1.0/RAND_MAX)*(kn_ll-2*topo_le)-(kn_ll-2*topo_le)/2)),ll));
      //
      //   cros.resize(1);
      //   cros[0].resize(2);
      //   displacements.resize(1)
      //   cros[0][0]=modmod((cns[0]-kn_ll/2.0),ll);
      //   cros[0][1]=modmod((cns[0]+kn_ll/2.0),ll);
      //
      //   type_sim.resize(1);
      //   type_sim[0]=1;
      //   kn_nn+=1;
      //
      //   }
      // }
      if(start==1){

      // cout<<displacements[0]<< " "<< displacements[1]<<endl;
      dx1[step]=displacements[0];
      dx2[step]=displacements[1];
    }
    }
    ////////////////
    ///END OF ONE SIM///
    ////////////////

    // for(int i=0;i<Tmax;i++) mean1+=dx1[i];
    // mean1=mean1/Tmax;
    // for(int i=0;i<Tmax;i++) mean2+=dx2[i];
    // mean2=mean2/Tmax;
    // double var1=0;
    // double var2=0;
    // for(int i=0;i<Tmax;i++) var1+=pow(dx1[i]-mean1,2);
    // for(int i=0;i<Tmax;i++) var2+=pow(dx2[i]-mean2,2);
    //
    // double den=pow(var1*var2,0.5);
    // for(int i=0;i<Tmax;i++) cor_coef+=(dx1[i]-mean1)*(dx2[i]-mean2)/den;
    // cout<<"Pearson correlation coefficient "<<cor_coef<<endl;
    //
    // for(int k=0;k<Tmax;k++) cor[k]+=dx1[k]*dx2[k];

    for(int k=0;k<Tmax;k++) {boo[k]+=knotted[k]; boolle[k][j]+=knotted[k];}



  }

  //Save data for post-processing
  for(int k=0;k<Tmax;k++) prob[k]=boo[k]*1.0/sims;

  // cout<<aaa<<endl;
  ///////////////////
  ///Writing files///
  ///////////////////

  // stringstream booleans;
  // booleans<<"booleans.dat";
  // ofstream writeboo;
  // writeboo.open(booleans.str().c_str());
  // for(int k=0;k<Tmax;k++) {
  //   for(int j=0;j<sims;j++) writeboo << boolle[k][j] << " ";
  //   writeboo << endl;
  // }
  // writeboo.close();

  cout<<"END"<<endl;

  //Export distances between knots to test free energy results
  // stringstream dista;
  // dista<<"distances.dat";
  // ofstream writedis;
  // writedis.open(dista.str().c_str());
  // for(int l=0;l<(sims);l++){
  // for(int k=0;k<Tmax;k++) {
  //   writedis << distances_sims[l][k] << endl;
  // }
  // }
  // writedis.close();

  cout<<"END"<<endl;
  ostringstream ss;
  if (bb==0){
    ss << "diff";
  }
  else{
    ss << "inter";
  }

  string title="pos_"+ ss.str()+".dat";

  stringstream posit;
  posit<<title;
  ofstream writeposit;
  writeposit.open(posit.str().c_str());

  for(int k=Tmax/2-5000;k<Tmax/2+5000;k++) {
    writeposit <<time_arr[k]<<" "<< pos1[k] <<" "<< endl;
  }

  writeposit.close();


  // for(int k=0;k<Tmax;k++) cor[k]=cor[k]/sims;
  // string cor_tit="correlation_"+ ss.str()+".dat";
  // stringstream cor_fil;
  // cor_fil<<cor_tit;
  // ofstream writecor;
  // writecor.open(cor_fil.str().c_str());
  //
  // for(int k=0;k<Tmax;k++) {
  //   writecor <<time_arr[k]<<" "<< cor[k]<< endl;
  // }
  // writecor.close();




  // string title="stds_"+ ss.str()+".dat";
  //
  // stringstream pos_stds;
  // pos_stds<<title;
  // ofstream writestds;
  // writestds.open(pos_stds.str().c_str());
  //
  // for(int k=0;k<Tmax;k++) {
  //   writestds <<time_arr[k]<<" "<< stds[k]<< endl;
  // }
  // writestds.close();


  // ostringstream ss;
  // ss << kn_nn_ori;
  // string aaa="prob_vs_time1_kn"+ ss.str()+".dat";
  //
  // stringstream prob_vs_time;
  // prob_vs_time << aaa;
  // ofstream writeprob;
  // writeprob.open(prob_vs_time.str().c_str());
  //
  // for(int k=0;k<Tmax;k++) writeprob << time_arr[k] << " " << prob[k] <<endl;
  // writeprob.close();

return 0;
}


/////////////////
///END OF MAIN///
/////////////////



double modmod(double num,double ln){

  return fmod(fmod(num,ln)+ln,ln);
}

double probab(double (&topo_pr)[knots_cons][knots_cons],double ran,int typ){
  double x=0.0;
  int j=0;
  for(int i=0;i<knots_cons;i++){
    x+=topo_pr[i][typ];
    if (x>=ran) break;
    j+=1;

  }
  return j;
}

double theta(double x){
  if(x>0) return 1;
  if(x<0) return 0;
}
double ff(double x){
  return (aa*x*x*(cos(pi*x/lzer)+1.0));
}
double fr(double x){

  return(bb*(ff(x)*(theta(lb-x)+theta(x-rb))+theta(x-lb)*theta(rb-x)*(ff(lb)+m*x-m*lb)));
}
double stand(double x[],int cut) {
  double mean=0;
  double stad=0;
  for(int i=0;i<cut;i++) mean+=x[i];
  mean=mean/cut;
  for(int i=0;i<cut;i++) stad+=pow((x[i]-mean),2);
  stad=sqrt(stad/cut);
  return stad;

}
