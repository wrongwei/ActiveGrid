/*
 Library to create sets of angles and angular speeds, in order to give to the grid a specific motion.
 
 There are a lot of functions, but also a lot of possible motions.
 
 The functions setanglesII and setspeeds are called from activegrid, and transform the sets of values (computed in algo) into orders for the grid.
 
 Florent Lachaussée
 March 2013
 
 Temporal correlation feature (aka 3D correlation) added by Kevin Griffin and Nathan Wei, summer 2015
 
 Dependencies: activegrid.cpp/.h, pickCorrelations.cpp/.h, loaf.cpp/.h
 
 */


#include <algo.h>
#include <assert.h>
#include <time.h>
#include <curses.h>

// modulo function with NON NEGATIVE remainder
inline int algo::modulo(int numerator, int denominator){
    int remainder;
    if (numerator*denominator<0){
        remainder = numerator%denominator - denominator*(floor(numerator/denominator)-1);
    }
    else {remainder=numerator%denominator;}
    
    return remainder;
}

//check the position of the servo and eventually switch angle from +angle to -angle
double  algo::checkorientationofservo(int col, int row, double angle){
    if((row%2==0 && col%2==0) || row%2==0)angle=-angle;
    return angle;
}


// as said in the name of the function
int  algo::setsameangletoallservos(double angle){
    double newangle[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0) newangle[col][row]= /*2*(floor(2*(float)rand()/RAND_MAX)- 0.5)**/checkorientationofservo(col,row,angle);
            else newangle[col][row]= Fictive;
        }
    }
    grid.setanglesII(newangle);
    return 1;
}



// clear name, used in allperiodic
int  algo::setsamespeedtoallservos(double angleperstep){
    //cout << "error1" << endl;
    double newangleperstep[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0){
                newangleperstep[col][row]=angleperstep;
            }
            else{
                newangleperstep[col][row]=0;
            }
        }
    }
    //cout << "error2" << endl;
    grid.setspeeds(newangleperstep);
    return 1;
}


//the name is clear
int  algo::setangleofoneservo(int tunedcol, int tunedrow, double tunedangle){
    double newangle[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0){
                if(col==tunedcol && row==tunedrow)newangle[col][row]=tunedangle;
                else newangle[col][row]=grid.angle[col][row];
            }
            else newangle[col][row]=Fictive;
        }
    }
    grid.setanglesII(newangle);
    return 1;
}


//name is clear
int  algo::setanglesofonerow(int tunedrow, double tunedangle){
    double newangle[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0){
                if(row==tunedrow)newangle[col][row]=tunedangle;
                else newangle[col][row]=grid.angle[col][row];
            }
            else newangle[col][row]=Fictive;
        }
    }
    grid.setanglesII(newangle);
    return 1;
}

//name is clear
int  algo::setanglesofonecolumn(int tunedcol, double tunedangle){
    double newangle[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0){
                if(col==tunedcol)newangle[col][row]=tunedangle;
                else newangle[col][row]=grid.angle[col][row];
            }
            else newangle[col][row]=Fictive;
        }
    }
    grid.setanglesII(newangle);
    return 1;
}

//closes the bottom-half, opens the top-half
int  algo::openhalfclosehalf(){
    double newangle[14][12];
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            if(grid.servo[col][row]!=0){
                if(row<=5)newangle[col][row]=90;
                else newangle[col][row]=0;
            }
            else newangle[col][row]=Fictive;
        }
    }
    grid.setanglesII(newangle);
    return 1;
}

// to set differents angles to all servos (several combinations between neighbours are possible)
int  algo::setanglestoallservos(float * positions, float * anglesteps, int combine, int constant, float rms){
    double newangle[14][12];
    double newangleperstep[14][12];
    int count=0;
    for(int row=1;row<12;row++){
        for(int col=1;col<14;col++){
            if(grid.servo[col][row]!=0){
                newangle[col][row]=positions[count];
                newangleperstep[col][row]=anglesteps[count];
                if(combine==1) count++;
                if(combine==2){
                    if (row==1 || col==13) count++;
                    if (row%2==1 && row!=1 && col!=13) {
                        if (col%2==1) {
                            newangle[col][row]=newangle[col+1][row-1];
                            newangleperstep[col][row]=newangleperstep[col+1][row-1];
                            count++;
                        }
                        if (col==12) count++;
                    }
                    if (row%2==0 && col%2==1) {
                        if (col!=1) {
                            newangle[col][row]=newangle[col-1][row-1];
                            newangleperstep[col][row]=newangleperstep[col-1][row-1];
                            count++;
                        }
                        else count++;
                    }
                }
                
                if (combine==3) {
                    if (row%4==1) {
                        if (col%3==0) {
                            newangle[col][row+1]=positions[count];
                            newangleperstep[col][row+1]=anglesteps[count];
                        }
                        if (col%3==1) {
                            if (row==5 && col==1) {
                                newangle[col][row-1]=positions[count];
                                newangleperstep[col][row-1]=anglesteps[count];
                            }
                            count++;
                        }
                        if (col%3==2) {
                            newangle[col-1][row+1]=positions[count];
                            newangleperstep[col-1][row+1]=anglesteps[count];
                            newangle[col][row+1]=positions[count];
                            newangleperstep[col][row+1]=anglesteps[count];
                            count++;
                        }
                        if (col==1) {
                            newangle[col][row-1]=positions[count];
                            newangleperstep[col][row-1]=anglesteps[count];
                            count++;
                        }
                        if (col==13) row++;
                    }
                    if (row%4==3) {
                        if (col%3==0) {
                            if (row!=11) {
                                newangle[col][row+1]=positions[count];
                                newangleperstep[col][row+1]=anglesteps[count];
                                newangle[col+1][row+1]=positions[count];
                                newangleperstep[col+1][row+1]=anglesteps[count];
                            }
                            count++;
                        }
                        if (col%3==2) {
                            if (row!=11) {
                                newangle[col][row+1]=positions[count];
                                newangleperstep[col][row+1]=anglesteps[count];
                            }
                            count++;
                        }
                        if (col==13) {
                            newangle[col][row-1]=positions[count];
                            newangleperstep[col][row-1]=anglesteps[count];
                            count++;
                        }
                        if (col==13) row++;
                    }
                }
                
                if(combine==4 || combine==8){
                    
                    newangle[col][row+1]=positions[count];
                    if(combine==8) newangle[col][row+1]=-positions[count];
                    newangleperstep[col][row+1]=anglesteps[count];
                    
                    if((col+(row%4))%2==1)count++;
                    if(col==13) {
                        row++;
                    }
                }
                if(combine==5){
                    if(col%5==1)count++;
                    if(col==12 || col==13 || col==1 || row==11 || (col==3 && row==10) || (col==4 && row==10) || (col==5 && row==10) || (col==6 && row==10)){
                        newangle[col][row]=positions[count];
                        newangleperstep[col][row]=anglesteps[count];
                    }
                }
                if(combine==6){
                    if(col%6==1)count++;
                    if(col==1 || col==13 || col==1 || row==11 || (col==8 && row==1) || (col==9 && row==1) || (col==10 && row==1) || (col==11 && row==1) || (col==12 && row==1) || (col==8 && row==9) || (col==9 && row==9) || (col==10 && row==9) || (col==11 && row==9) || (col==12 && row==9) || row==10 || row==11){
                        newangle[col][row]=positions[count];
                        newangleperstep[col][row]=anglesteps[count];
                    }
                }
                if(combine==9){
                    if((col+(row%4))%2==1) {
                        newangle[col-1][row+1]=positions[count];
                        newangleperstep[col-1][row+1]=anglesteps[count];
                        count++;
                    }
                    if((col+(row%4))%2==0){
                        newangle[col][row]=-positions[count];
                        if(col!=13) {
                            newangle[col+1][row+1]=-positions[count];
                            newangleperstep[col+1][row+1]=anglesteps[count];
                        }
                    }
                    if(col==13) {
                        newangle[col][row+1] = -positions[count];
                        newangleperstep[col][row+1] = anglesteps[count];
                        row++;
                    }
                }
            }
            
            
            else{ // for non-existing servos
                if(combine==4 || combine==8 || combine==9){
                    if (row==1&&col==1) {
                        if(combine==4) newangle[col][row+1]=positions[count];
                        if(combine==8) newangle[col][row+1]=-positions[count];
                        newangleperstep[col][row+1]=anglesteps[count];
                        if(combine==9){
                            newangle[col+1][row+1]=positions[count];
                            newangleperstep[col+1][row+1]=anglesteps[count];
                        }
                    }
                    if (row==1&&col==13) {
                        newangle[col][row]=Fictive;
                        newangleperstep[col][row]=0;
                        newangle[col][row+1]=positions[count];
                        newangleperstep[col][row+1]=anglesteps[count];
                        count++;
                    }
                    if((col+(row%4))%2==1)count++;
                    if(col==13) {
                        row++;
                    }
                }
                if (combine==3) {
                    if ((row==1 && col==13) || (row==9 && col==13)) {
                        newangle[col][row]=Fictive;
                        newangleperstep[col][row]=0;
                        row++;
                        count++;
                    }
                    if (row==9 && col==1) {
                        newangle[col][row-1]=positions[count];
                        newangleperstep[col][row-1]=anglesteps[count];
                        count++;
                    }
                }
                if (row!=2) {
                    newangle[col][row]=Fictive;
                    newangleperstep[col][row]=0;
                }
            }
        }
    }
    if(constant==1) {area(newangle, rms);}
    
    grid.setspeeds(newangleperstep);
    grid.setanglesII(newangle);
    
    // writing-on-file can be commented to save computational time
    // for plot-output-file
    /*for(int row=1;row<12;row++){
     for(int col=1;col<14;col++){
     if(grid.servo[col][row]!=0){
     anglefile << "    " << newangle[col][row];}
     }
     }
     anglefile << endl;*/
    
    return 1;
}


// to set differents angles to all servos (without combination between neighbours)
// light version for correlated motion (less computational time), used in all spatial correlation methods
int algo::setanglestoallservosII(float * positions, float * anglesteps, int constant, float rms){
    double newangle[14][12];
    double newangleperstep[14][12];
    int count=0;
    for(int row = 1; row < 12; row++) {
        for(int col = 1; col < 14; col++) {
            if (grid.servo[col][row]!=0){
                newangle[col][row]=positions[count];
                newangleperstep[col][row]=anglesteps[count];
                count++;
            }
            else { // for non-existing servos
                if (row!=2) {
                    newangle[col][row]=Fictive;
                    newangleperstep[col][row]=0;
                }
            }
        }
    }
    
    if(constant==1) {area(newangle, rms);}
    
    grid.setspeeds(newangleperstep);
    grid.setanglesII(newangle);
    
    // writing-on-file removed for space; see above method for angle-writing code
    
    // write all angles to file
    for(int i = 0; i < 143; i++){
        anglefile << "    " << positions[i];
    }
    anglefile << endl;
    
    return 1;
}

// give the same periodic motion to all paddles
int algo::allperiodic(double angle, double frequency){
    
    //timing:
    timeval testtime;
    gettimeofday(&testtime,0);
    cout << "sec: " << testtime.tv_sec << endl;
    cout << "microsec: " << testtime.tv_usec << endl;
    
    int numberofloops=20;
    int steps=(int)((1000000/(double)updatetimeinmus)*frequency);
    double angleperstep=(angle*4/steps);
    cout << "frequency: " << frequency << " [Hz], updaterate: " << (1000000/(double)updatetimeinmus) << " [Hz], steps for half period: " << steps/2 << "\nangle: " << angle << ", angleperstep: " << angleperstep << endl;
    
    double actualangle=angle*(-1);
    int sign=-1;
    setsameangletoallservos(actualangle);
    
    gettimeofday(&testtime,0);
    long time_usec=0;//testtime.tv_usec;
    //  long time_sec=testtime.tv_sec;
    while ( testtime.tv_usec > updatetimeinmus ) gettimeofday(&testtime,0);
    
    
    //loop over number of periodic movements:
    for(int count = 0;count<=(numberofloops*2);count++){
        
        cout << "\nnumber of loops: " << count << endl;
        cout << "actual sec: " << testtime.tv_sec << endl;
        cout << " angle: " << endl;
        
        for(int anglestep = 1;anglestep<=(steps/2);anglestep++){
            time_usec += updatetimeinmus;
            gettimeofday(&testtime,0);
            if(time_usec>1000000)time_usec-=1000000;
            cout << "difference: " << (time_usec - testtime.tv_usec)/1000 << " milli sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
            
            if(testtime.tv_usec > time_usec) {
                cout << "---------------------------------Problem!!!------------------------------" << endl;
                //cout << "difference: " << (time_usec - testtime.tv_usec) << " mu sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
            }
            
            while (testtime.tv_usec <= time_usec){
                gettimeofday(&testtime,0);
                if(time_usec==1000000 && testtime.tv_usec<updatetimeinmus ){
                    break;
                }
            }
            
            cout << actualangle << ", " << endl;
            
            if(actualangle>=(angle)){
                sign=-1;
            }
            if(actualangle<=(angle*(-1))){
                sign=+1;
            }
            
            actualangle+=angleperstep*sign;
            
            setsameangletoallservos(actualangle);
            setsamespeedtoallservos(angleperstep);
            
        }
    }
    
    
    return 0;
}


// give a chaotic movement to the paddles. Combinations between paddles are possible
int algo::chaoticMovement(int combine, int constant, int option){
    
    float rms =0;
    
    anglefile.open("angleservo2.txt", ios::out | ios::trunc); // file to plot angles in function of time
    for (int numero=0; numero < 129; numero++){
        anglefile << "   Angle(" << numero << ")";}
    anglefile << endl;
    
    // takes a first random correlated sequence of angles with the same parameters
    // compute its rms value of angles, which will be used to keep the area constant
    rms=compute_rms(option);
    
    // initialize values to 0 (necessary to avoid NaN output)
    for (int i=0; i<numberOfServos; i++){
        positions[i]=0;
        anglesteps[i]=0;
    }
    
    //timing:
    timeval testtime;
    gettimeofday(&testtime,0);
    long time_usec=0;
    while ( testtime.tv_usec > updatetimeinmus) gettimeofday(&testtime,0);
    
    
    // main loop : regularly calls "run" to have new positions and angles, and send them to "setanglestoallservos"
    while(0==0){
        
        //getpositions of each servo:
        run(positions,anglesteps,option);
        
        //setposition of each servo:
        time_usec += updatetimeinmus;
        gettimeofday(&testtime,0);
        if(time_usec>1000000)time_usec-=1000000;
        if(testtime.tv_usec > time_usec) {
            cout << "---------------------------------Problem!!!------------------------------" << endl;
            cout << "difference: " << (time_usec - testtime.tv_usec) << " mu sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
        }
        while (testtime.tv_usec <= time_usec){
            gettimeofday(&testtime,0);
            if(time_usec==1000000 && testtime.tv_usec<updatetimeinmus ){
                break;
            }
        }
        setanglestoallservos(positions,anglesteps, combine ,constant, rms); // for motion
        
    }
    anglefile.close();
    return 0;
}



// Takes a first random correlated sequence of angles with the same parameters as in input
// compute its rms value of angles, which will be used to keep the area constant.
float algo::compute_rms(int option){
    float control_positions[numberOfServos][4000];
    float mean=0;
    float rms =0;
    
    // first loop to compute the rms value of angles, which will be used to keep the area constant
    // takes a first random correlated sequence of angles, without correction
    for (int t =0; t<4000;t++){
        run(positions,anglesteps,option);
        
        for (int i=0; i<numberOfServos;i++){
            control_positions[i][t]=positions[i];
        }
    }
    
    // computes mean and rms value of angles in the first sequence
    for (int i=0; i<numberOfServos;i++){
        for (int t=0;t<4000;t++){
            mean+= control_positions[i][t]/(numberOfServos*4000);
        }
    }
    for (int i=0; i<numberOfServos;i++){
        for (int t=0;t<4000;t++){
            rms+= pow(control_positions[i][t]-mean, (int)2)/(numberOfServos*4000);
        }
    }
    
    rms=sqrt (rms); // rms is the sqrt of variance (that we actually computed)
    
    return rms;
}


// reads the contents of positions and steps, put the value stored at actualposition in vector in actpos and actstep, and actualpositioninvector++
void algo::run(float actpos[], float actstep[], int option){
    for(int i=0;i<numberOfServos;i++){
        //get current positions:
        if(positions_random[i].at(actualpositioninvector[i])!=-100){
            actpos[i]=positions_random[i].at(actualpositioninvector[i]);
            actstep[i]=steps_random[i].at(actualpositioninvector[i]);
            actualpositioninvector[i]++;
        }
        else{ // compute the next steps
            if (option==1) updateOneWing(i); // for piecewise periodic motion (choice 11)
            else updateOneWing2(i); // for random motion (choice 12)
            actpos[i]=positions_random[i].at(actualpositioninvector[i]);
            actstep[i]=steps_random[i].at(actualpositioninvector[i]);
            actualpositioninvector[i]++;
        }
    }
}




// generates a piecewise periodic random sequence of orders for one paddle and store them in positions and steps.
void algo::updateOneWing(int WingNumber){
    //clear positions:
    positions_random[WingNumber].clear();
    steps_random[WingNumber].clear();
    
    //calculate randomly speed and maximum amplitude for next step
    speed=(((float)rand()/RAND_MAX)*(max_speed-min_speed))+min_speed;
    amplitude=(((float)rand()/RAND_MAX)*(max_amplitude-min_amplitude))+min_amplitude;
    
    //it might occur that the calculated speed is higher than the maximum amplitude. In this case the speed is adapted to the amplitude to avoid errors.
    if (speed > amplitude) speed=amplitude;
    
    //calculate frequency for user:
    frequency=5*speed/(2*amplitude);
    //cout << frequency << endl;
    
    //calculate posisitons:
    double actualangle=0;
    int sign=+1;
    
    //entire sequence for one servo is calculated including positions and steps
    for(int k=0;k<(int)(floor((4*amplitude/speed)+0.5));k++){ //Total number of steps
        if((actualangle+speed)>=(amplitude)){
            sign=-1;
            actualangle=amplitude - (speed-(amplitude-actualangle));
        }
        if((actualangle-speed)<=(amplitude*(-1))){
            actualangle=(-1)*(amplitude-(speed-(amplitude+actualangle)));
            sign=+1;
        }
        if(actualangle>(amplitude*(-1)) && actualangle<(amplitude)){
            actualangle+=speed*sign;
        }
        //cout << actualangle << ", ";
        if(positions_random[WingNumber].size()>0){
            if(positions_random[WingNumber].at(positions_random[WingNumber].size()-1)<0 && actualangle>0) break;
            else {
                positions_random[WingNumber].push_back(actualangle);
                steps_random[WingNumber].push_back(speed);
            }
        }
        else{
            positions_random[WingNumber].push_back(actualangle);
            steps_random[WingNumber].push_back(speed);
        }
    }
    //cout << endl;
    //end:
    positions_random[WingNumber].push_back(-100);
    steps_random[WingNumber].push_back(-100);
    actualpositioninvector[WingNumber]=0;
}


// generates a "totally" random sequence of angles for one paddle and store orders (angles and speeds) in positions and steps.
void algo::updateOneWing2(int WingNumber){
    //clear positions:
    positions_random[WingNumber].clear();
    steps_random[WingNumber].clear();
    
    //calculate positions:
    double actualangle;
    double Amplitude;
    int nbofsteps;
    
    actualangle=old_angle[WingNumber];
    
    //calculate randomly angle for next step
    new_angle[WingNumber]=(((float)rand()/RAND_MAX)*(max_angle-min_angle))+min_angle;
    
    // calculate the amplitude between the old and the new angles
    Amplitude = new_angle[WingNumber] - old_angle[WingNumber];
    
    // calculate the number of steps necessary to go through the amplitude range without exceeding maximal speed
    nbofsteps = (1 + floor(fabs(Amplitude)/max_speed));
    
    // compute speed
    speed = Amplitude/(nbofsteps);
    
    //entire sequence for one servo is calculated including positions and steps
    for(int k=0;k< nbofsteps ;k++){ //Total number of steps
        actualangle += speed;
        // put nbofsteps entries onto positions and steps vectors
        positions_random[WingNumber].push_back(actualangle);
        steps_random[WingNumber].push_back(speed);
    }
    //cout << endl;
    //end:
    positions_random[WingNumber].push_back(-100); // mark end of sequence with -100 (so system knows angle has been reached)
    steps_random[WingNumber].push_back(-100);
    actualpositioninvector[WingNumber]=0;
    old_angle[WingNumber]=new_angle[WingNumber]; // == actualangle
}


// chaotic motion with correlation between paddles
int algo::correlatedMovement(int constant, float sigma, float alpha, double height, int mode, int mrow, int mcol, float target_rms){
    
    float rms;
    float correction=1;
    int i=0;
    
    // file header for just the angles that get sent to actual servos
    /* anglefile.open("angleservo2.txt", ios::out | ios::trunc); // file to plot angles in function of time
     for (int numero=0; numero < 129; numero++){
     anglefile << "   Angle(" << numero << ")";
     }
     anglefile << endl; */
    
    // file to plot angles in function of time
    // UI to avoid accidentally overwriting angle files
    ifstream ifile("angleservo_cM.txt");
    int overwriteFile;
    if (ifile) {
        cout << "WARNING: An angle file for this program already exists. 1: continue and overwrite file, 0: kill program" << endl;
        cin >> overwriteFile;
        if (!overwriteFile)
            exit(0);
    }
    anglefile.open("angleservo_cM.txt", ios::out | ios::trunc);
    /*for (int numero=0; numero < 143; numero++){
     anglefile << "   Angle(" << numero << ")";
     }
     anglefile << endl;*/
    
    /* Declare function pointers for the spatial correlation functions
     First, we declare a pointer *pfCorr to a function with 2 arguments j and k, plus associated parameters needed for correlation function operation (currently just sigma and height).
     Then we call the function pickSpatialCorr (from pickCorrelations.cpp).
     This function returns a pointer to the function of choice (determined by 'mode'),
     which can then be used by simply calling pfCorr(j, k). */
    float (*pfCorr)(int j, int k, float sigma, float alpha, float height);
    pfCorr = pickSpatialCorr(mode);
    
    // compute normalization for any convolution formula
    norm1 = 0; // stored in algo.h (for accessibility)
    float bound; // set bounds of correlation kernel
    if (mode <= 4) bound = range_of_corr; // = 7, defined in algo.h
    else bound = sigma;
    // Note: this is different from the previous implementation, which had mysteriously different logic for each function
    if (mode == 6 || mode == 7) norm1 = 1; // no normalization necessary for modes 6 and 7
    else {
        for (int j = -bound; j <= bound; j++) { // range of neighbors used to compute normalization/convolution
            for (int k = -bound; k <= bound; k++) { // j and k refer to the shift
                norm1 += pfCorr(j, k, sigma, alpha, height); // normalization is sum of values of kernel function across the range of the kernel
            }
        }
    }
    
    // makes a random correlated sequence of angles, with the same parameters but without correction
    // computes its mean and rms value of angles. This is done so that the rms correction factor can be
    // determined before the angles have been produced
    rms=compute_rmscorr(sigma, mode, alpha, height, mrow, mcol);
    
    
    correction=target_rms/rms; // correction factor
    
    // initialize values to 0 (necessary to avoid NaN output)
    for (int i=0; i<numberOfServos; i++){
        positions[i]=0;
        anglesteps[i]=0;
        old_positions[i]=0;
        old_steps[i]=0;
        err[i]=0;
    }
    
    //timing:
    timeval testtime;
    gettimeofday(&testtime,0); // writes the time into the feilds of the structure named testtime
    long time_usec=0; // set a timer named time_usec to zero
    while ( testtime.tv_usec > updatetimeinmus) gettimeofday(&testtime,0); // keep updating the usec field of the absolute time, while the testtime is > 0.1 seconds. So at most after .9 seconds exit the loop
    
    // main loop: give angle orders
    while(0==0){
        
        //getpositions of each servo: (i.e. run correlation functions using runcorr method)
        runcorr(positions,anglesteps,sigma,alpha,height,mode,mrow,mcol,correction,norm1,old_positions,old_steps,err);
        
        //setposition of each servo:
        time_usec += updatetimeinmus;  // add the 0.1 sec to the stopwatch
        gettimeofday(&testtime,0);   // update the abolute time
        if(time_usec>1000000)time_usec-=1000000; // if the stopwatch is greater than 1 sec, subtract 1 sec
        if(testtime.tv_usec > time_usec) { // if the absolute time > the stopwatch, there is a problem (runcorr took too long)
            cout << "---------------------------------Problem!!!------------------------------" << endl;
            cout << "difference: " << (time_usec - testtime.tv_usec) << " mu sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
        }
        while (testtime.tv_usec <= time_usec){// while the absolute time <= the stopwatch
            gettimeofday(&testtime,0); // update the absolute time
            if(time_usec==1000000 && testtime.tv_usec<updatetimeinmus ){ // if the stopwatch equals 1 sec && the absolute time < 0.1 sec, break
                break;
            }
        }
        setanglestoallservosII(positions,anglesteps,constant,target_rms); // for motion - send resulting angles stored in positions vector to the grid
        cout << i << "\n"; // print index of iteration (for UI purposes)
        i += 1;
        
    }
    anglefile.close();
    return 0;
}


// take a random sequence and computes its std dev. It's useful for the correction coefficent that
// is needed to give to the output the desired rms value of angles.
float algo::compute_rmscorr(float sigma, int mode, float alpha, double height, int mrow, int mcol){
    
    float control_positions[numberOfServos][4000];
    float mean=0;
    float rms=0;
    
    // takes a random correlated sequence of angles, without correction
    for (int t =0; t<4000;t++){
        runcorr(positions,anglesteps,sigma,alpha,height,mode,mrow,mcol,1,norm1,old_positions, old_steps, err);
        
        for (int i=0; i<numberOfServos;i++){
            control_positions[i][t]=positions[i];
        }
    }
    
    // computes mean and rms value of angles in the first sequence
    for (int i=0; i<numberOfServos;i++){
        for (int t=0;t<4000;t++){
            mean+= control_positions[i][t]/(numberOfServos*4000);
        }
    }
    for (int i=0; i<numberOfServos;i++){
        for (int t=0;t<4000;t++){
            rms+= pow(control_positions[i][t]-mean, (int)2)/(numberOfServos*4000);
        }
    }
    
    rms=sqrt (rms); // rms is the sqrt of variance (that we actually computed)
    return rms;
    
}

// reads the contents of positions and steps, put the value stored at actualposition in vector in actpos and actstep, and actualpositioninvector++
// then creates correlation between paddles
void algo::runcorr(float actpos[], float actstep[], float sigma, float alpha, double height, int mode, int mrow, int mcol, float correction, float norm, float oldpos[], float oldstep[], float err[]){
    
    float interpos[13][11]; // arrays storing the intermediate computed values before convolution
    //float interstep[13][11]; // unused
    
    for(int i=0;i<numberOfServos;i++){
        //get current positions:
        if(positions_random[i].at(actualpositioninvector[i])==-100){ // if the end of the positions vector is reached
            updateOneWing2(i); // get more random angles for this servo from updateOneWing2, along with steps
        }
        interpos[columns[i]][rows[i]]=positions_random[i].at(actualpositioninvector[i]); // convert these angles from vector form into 2d array form
        //interstep[columns[i]][rows[i]]=steps_random[i].at(actualpositioninvector[i]); // unused
        actualpositioninvector[i]++; // move position counter up one, to mark next angle (for next iteration)
    }
    // Declare function pointers for the spatial correlation functions (see note in correlatedMovement)
    float (*pfCorr)(int j, int k, float sigma, float alpha, float height);
    pfCorr = pickSpatialCorr(mode);
    
    // convolution to create correlation between paddles
    // periodic boundary conditions are used
    
    // preliminary computations for special cases
    float bound; // determine bounds of iteration for for loops
    bool isGaussian = false;
    if (mode == 1) isGaussian = true;
    if (mode <= 4) bound = range_of_corr;
    else bound = sigma;
    int col = 0; int row = 0;
    if (mode == 7) {
        int r = rand() % 143; // generate random paddle
        mcol = modulo(columns[r],13); // set mcol and mrow to denote random paddle
        mrow = modulo(rows[r],11);
    }
    
    // Loop through servos and calculate/create correlations, using helper methods (below)
    // for a given servo: take the dot product of the kernel (centered at the given paddle) with
    // the angles lying within the range of the kernel, then normalize the result to get new correlated angle
    for (int i = 0; i < numberOfServos; i++) {
        
        initialize_pos_step(actpos, actstep, oldpos, oldstep, i); // set up old/act arrays
        
        // true top hat (not modular yet)
        if (mode == 6 || mode == 7) {
            actpos[i] = 0; // keep it functional for now
            truetophat2d(actpos, oldpos, actstep, correction, interpos, sigma, mrow, mcol, i);
        }
        // master loop for modular functions (can remove else if true top hat becomes modularized)
        else {
            for (int j = -bound; j <= bound; j++) { // range of neighbours used to compute convolution
                for (int k = -bound; k <= bound; k++) { // j and k refer to the shift
                    col = modulo(columns[i] + j, 13); // controls periodic boundary conditions
                    row = modulo(rows[i] + k, 11);
                    actpos[i] += correction * pfCorr(j, k, sigma, alpha, height) * interpos[col][row];
                }
            }
            actpos[i] = actpos[i] / norm; // normalize by pre-calculated coefficient
        }
        
        if (mode == 1) actstep[i]= actpos[i]-oldpos[i]+err[i];
        else actstep[i] = actpos[i]-oldpos[i];
        
        safety_check(actpos, actstep, err, isGaussian, i); // angle-checking function
    }
}

// Correlation helper methods

// initializes arrays to zero
void algo::initialize_pos_step(float actpos[], float actstep[], float oldpos[], float oldstep[], int i) {
    oldpos[i]=actpos[i];
    oldstep[i]=actstep[i]; // oldstep is currently not used
    actpos[i]=0; // clear the output
    actstep[i]=0;
}

// modifies arrays to reflect true top hat correlation (random or user-specified) in 2d
void algo::truetophat2d(float actpos[], float oldpos[], float actstep[], float correction,
                        float interpos[][11], float sigma, int mrow, int mcol, int i) {
    int col = modulo(columns[i],13);
    int row = modulo(rows[i],11);
    
    double dist = sqrt(((col-mcol)*(col-mcol))+((row-mrow)*(row-mrow)));
    if (dist <= sigma){
        actpos[i] = correction * interpos[mcol][mrow];
    } else {
        actpos[i] = correction * interpos[col][row];
    }
}

// Makes sure angles and speeds are not out of the range of servo motor capabilities
void algo::safety_check(float actpos[], float actstep[], float err[], bool isGaussian, int i) {
    // angle safety: do not exceed amplitude of 90 degrees
    if (actpos[i]>90) {
        actpos[i]=90;
        actstep[i]=0;
    }
    else if (actpos[i]<-90) {
        actpos[i]=-90;
        actstep[i]=0;
    }
    // speed safety (Gaussian): do not exceed servo speed of 60 degrees in 0.15 seconds
    if (isGaussian) {
        if (actstep[i] > max_speed) {
            err[i]=actstep[i] - max_speed;
            actstep[i] = max_speed;
        }
        else if (actstep[i] < -max_speed) {
            err[i]=actstep[i] + max_speed;
            actstep[i] = -max_speed;
        }
    }
    // speed safety (non-Gaussian): same as above, without err array
    else {
        if (actstep[i] > max_speed)
            actstep[i] = max_speed;
        else if (actstep[i] < -max_speed)
            actstep[i] = -max_speed;
    }
}

// function to keep the projected area of the grid constant
// (NOTE: this is known to cause grid freezing issues when running 3D correlations with constant area)
void algo::area(double actpos[14][12], float rms) {
    //second attempt to keep the projected area constant
    angle_constant_area = 129 * (100*sin((rms/1.22) * (M_PI / 180))); // 24.4779
    double projected_area=0;
    //T_1 and T_2 are coefficients that need to be calculated
    double T_1=0;
    double T_2=0;
    double theta_correct;
    for (int row=1; row<12; row++) {
        for (int col=1; col<14; col++) {
            if (actpos[col][row]!=Fictive) {
                T_1 += fabs(sin(actpos[col][row] * M_PI / 180));
                T_2 += fabs(cos(actpos[col][row] * M_PI / 180));
            }
        }
    }
    //cout << "T_1 :" << T_1 << "  T_2 :" << T_2 << endl;
    //cout << "constant area :" << angle_constant_area << endl;
    //theta_correct is the correction angle for all servos to keep the projected area constant
    theta_correct = (asin(angle_constant_area / 100 / sqrt(pow(T_1,2) + pow(T_2,2))) + atan(T_2/T_1) - M_PI/2)/M_PI*180;
    //cout << theta_correct << endl;
    //calculate the projected area before correction
    for (int row=1; row<12; row++) {
        for (int col=1; col<14; col++) {
            if (actpos[col][row]!=Fictive) {
                projected_area += 100 * fabs(sin(actpos[col][row] * M_PI / 180));
            }
            //cout << j << "  before:" << "actpos: " << actpos[j] << endl;
        }
    }
    //cout << "old projected area: " << projected_area << endl;
    //apply correction for current step
    for (int row=1; row<12; row++) {
        for (int col=1; col<14; col++) {
            if (actpos[col][row]!=Fictive) {
                if (actpos[col][row]<0) {
                    actpos[col][row]-=theta_correct;
                    if (actpos[col][row]>0) actpos[col][row]=0;
                }
                else {
                    actpos[col][row]+=theta_correct;
                    if (actpos[col][row]<0) actpos[col][row]=0;
                }
            }
        }
        
    }
    
    // check if there is no angle angle greater than 90 degrees, correct them if necessary
    for (int row=1; row<12; row++) {
        for (int col=1; col<14; col++) {
            if (actpos[col][row]!=Fictive) {
                if (actpos[col][row]>90) actpos[col][row]=90;
                else if (actpos[col][row]<-90) actpos[col][row]=-90;
            }
        }
    }
    
    
    
    //calculate the projected area after application of the correction
    projected_area=0;
    for (int row=1; row<12; row++) {
        for (int col=1; col<14; col++) {
            if (actpos[col][row]!=Fictive) {
                projected_area += 100 * fabs(sin(actpos[col][row] * M_PI / 180));
            }
            //cout << j << "  after:" << "actpos: " << actpos[j] << endl;
        }
    }
    //cout << "new projected area :" << projected_area << endl;
    /*there are some exceptions when the algorithm allows peaks, in those cases the second correction
     reduces the peaks by a manual reduction of the angle. the servos are chosen randomly and it stops
     when a limit is crossed. */
    if (projected_area > (angle_constant_area + 45)) {
        int random[143];
        srand(time(0));
        for (int i=0; i<143; i++) {
            random[i]=i;
        }
        for (int n=0; n<143; n++) {
            int r = n + (rand() % (143 - n));
            int temp = random[n]; random[n] = random[r]; random[r] = temp;
        }
        for (int k=0; k < 143; k++) {
            int row = ceil(random[k]/13);
            int col = random[k] - (row * 13);
            if (fabs(actpos[col][row]) < 45 && fabs(actpos[col][row]) > 5) {
                if (actpos[col][row] > 0) actpos[col][row] -= 3;
                else actpos[col][row] += 3;
            }
            projected_area=0;
            for (int row=1; row<12; row++) {
                for (int col=1; col<14; col++) {
                    if (actpos[col][row]!=Fictive) {
                        projected_area += 100 * fabs(sin(actpos[col][row] * M_PI / 180));
                    }
                }
            }
            if (projected_area < (angle_constant_area + 20)) break;
        }
    }
}


// 2-step correlated motion between paddles (choice 14 in menuII)
int algo::correlatedMovement_steps(int constant, float sigma1, float sigma2, int mode, float target_rms1, float target_rms2, double period, double duty_cycle){
    
    float rms1, rms2;
    float correction1=1;
    float correction2=1;
    float target_rms; // used to keep the area constant
    double temps;
    long i=0;
    
    // UI to avoid accidentally overwriting angle files
    ifstream ifile("angleservo_cMs.txt");
    int overwriteFile;
    if (ifile) {
        cout << "WARNING: An angle file for this program already exists. 1: continue and overwrite file, 0: kill program" << endl;
        cin >> overwriteFile;
        if (!overwriteFile)
            exit(0);
    }
    anglefile.open("angleservo_cMs.txt", ios::out | ios::trunc); // file to plot angles in function of time
    /*for (int numero=0; numero < 129; numero++){
     anglefile << "   Angle(" << numero << ")";}
     anglefile << endl;*/
    
    // Declare function pointers for the spatial correlation functions
    float (*pfCorr)(int j, int k, float sigma, float alpha, float height);
    pfCorr = pickSpatialCorr(mode);
    
    // compute normalization for any convolution formula
    norm1=0;
    // Note: this is different from the previous implementation (but I'm not sure whether it works like the other one either)
    for (int j = -range_of_corr; j <= range_of_corr; j++) { // range of neighbors used to compute normalization/convolution
        for (int k = -range_of_corr; k <= range_of_corr; k++) { // j and k refer to the shift
            norm1 += pfCorr(j, k, sigma1, 0, 0);
        }
    }
    
    norm2=0;
    // Note: this is different from the previous implementation, which had mysteriously different logic for each function
    for (int j = -range_of_corr; j <= range_of_corr; j++) { // range of neighbors used to compute normalization/convolution
        for (int k = -range_of_corr; k <= range_of_corr; k++) { // j and k refer to the shift
            norm2 += pfCorr(j, k, sigma2, 0, 0);
        }
    }
    
    // takes a first random correlated sequence of angles, with the same parameters but without correction
    // computes its mean and rms value of angles
    rms1=compute_rmscorr(sigma1, mode, 0,0,0,0);
    correction1=target_rms1/rms1; // correction factor for first half-period
    
    rms2=compute_rmscorr(sigma2, mode, 0,0,0,0);
    correction2=target_rms2/rms2;
    
    // initialize values to 0 (necessary to avoid NaN output)
    for (int i=0; i<numberOfServos; i++){
        positions[i]=0;
        anglesteps[i]=0;
        old_positions[i]=0;
        old_steps[i]=0;
        err[i]=0;
    }
    
    //timing:
    timeval testtime;
    gettimeofday(&testtime,0);
    long time_usec=0;
    while ( testtime.tv_usec > updatetimeinmus) gettimeofday(&testtime,0);
    
    
    // main loop: give angle orders
    while(0==0){
        
        //getpositions of each servo:
        gettimeofday(&testtime,0);
        temps = testtime.tv_sec + (testtime.tv_usec/1000000.0);
        if ( ((temps/period)- floor(temps/period)) < duty_cycle ){
            runcorr(positions,anglesteps,sigma1,0,0,mode,0,0,correction1,norm1, old_positions, old_steps, err);
            target_rms=target_rms1;
            grid.high_duty = true;
        }
        else {
            runcorr(positions,anglesteps,sigma2,0,0,mode,0,0,correction2,norm2, old_positions, old_steps, err);
            target_rms=target_rms2;
            grid.high_duty = false;
        }
        
        //setposition of each servo:
        time_usec += updatetimeinmus;
        gettimeofday(&testtime,0);
        if(time_usec>1000000)time_usec-=1000000;
        if(testtime.tv_usec > time_usec) {
            cout << "---------------------------------Problem!!!------------------------------" << endl;
            cout << "difference: " << (time_usec - testtime.tv_usec) << " mu sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
        }
        while (testtime.tv_usec <= time_usec){
            gettimeofday(&testtime,0);
            if(time_usec==1000000 && testtime.tv_usec<updatetimeinmus ){
                break;
            }
        }
        setanglestoallservosII(positions,anglesteps,constant, target_rms); // for motion
        cout << i << "\n";
        i += 1;
        
    }
    
    return 0;
}

// chaotic motion with correlation between paddles (choice 15 in menuII)
int algo::correlatedMovement_periodic(int constant, float sigma, int mode, float target_rms, int numberofsteps){
    
    float rms;
    float correction=1;
    
    // tools to implement the periodicity follow
    int iteration=0;
    int phase = 0;
    int loop_number = 0;
    int j;
    float corrected_speed;
    bool smooth;
    float stored_positions[numberOfServos][numberofsteps];
    float stored_steps[numberOfServos][numberofsteps];
    
    // UI to avoid accidentally overwriting angle files
    ifstream ifile("angleservo_cMp.txt");
    int overwriteFile;
    if (ifile) {
        cout << "WARNING: An angle file for this program already exists. 1: continue and overwrite file, 0: kill program" << endl;
        cin >> overwriteFile;
        if (!overwriteFile)
            exit(0);
    }
    anglefile.open("angleservo_cMp.txt", ios::out | ios::trunc); // file to plot angles in function of time
    /*for (int numero=0; numero < 129; numero++){
     anglefile << "   Angle(" << numero << ")";}
     anglefile << endl;*/
    
    // Declare function pointers for the spatial correlation functions
    float (*pfCorr)(int j, int k, float sigma, float alpha, float height);
    pfCorr = pickSpatialCorr(mode);
    
    // compute normalization for any convolution function (necessary, even with the correction)
    norm1=0;
    for (int j = -range_of_corr; j <= range_of_corr; j++) { // range of neighbors used to compute normalization/convolution
        for (int k = -range_of_corr; k <= range_of_corr; k++) { // j and k refer to the shift
            norm1 += pfCorr(j, k, sigma, 0, 0);
        }
    }
    
    // takes a first random correlated sequence of angles, with the same parameters but without correction
    // computes its mean and rms value of angles
    rms=compute_rmscorr(sigma, mode, 0,0,0,0);
    
    correction=target_rms/rms; // correction factor
    
    // initialize values to 0 (necessary to avoid NaN output)
    for (int i=0; i<numberOfServos; i++){
        positions[i]=0;
        anglesteps[i]=0;
        old_positions[i]=0;
        old_steps[i]=0;
        err[i]=0;
    }
    
    
    // computes the random sequence of the periodic pattern
    for (int t =0; t<numberofsteps;t++){
        runcorr(positions, anglesteps, sigma, 0, 0, mode, 0, 0, correction, norm1, old_positions, old_steps, err);
        
        for (int i=0; i<numberOfServos;i++){
            stored_positions[i][t]=positions[i];
            stored_steps[i][t]=anglesteps[i];
        }
    }
    
    
    // corrects the end of the pattern to make a smooth connection with its beginning
    for (int i=0; i<numberOfServos;i++){
        j=1; // counts the number of final steps you have to modify
        smooth =false;
        while (smooth== false){
            corrected_speed = (stored_positions[i][0]-stored_positions[i][numberofsteps-j])/j; // speed necessary to reach the initial value in time
            if (corrected_speed <= 40.0){ // safety condition
                for (int k=1; k<j; k++){
                    stored_positions[i][numberofsteps-j+k] = stored_positions[i][numberofsteps-j]+ k * corrected_speed;
                    stored_steps[i][numberofsteps-j+k] = corrected_speed;
                }
                stored_steps[i][numberofsteps-j] = corrected_speed;
                smooth=true;
            }
            else {j++;}
        }
        
    }
    
    /*angleperiod.open("angleservo_periodic.txt", ios::out | ios::trunc); // file to plot angles in function of time
     for (int numero=0; numero < 129; numero++){
     angleperiod << "   Angle(" << numero << ")";}
     angleperiod << endl;
     
     for (int t=0; t<numberofsteps;t++){
     for (int i=0; i<129;i++){
     angleperiod << "    " << stored_positions[i][t];
     }
     angleperiod << endl;
     }
     
     angleperiod.close();*/
    
    //timing:
    timeval testtime;
    gettimeofday(&testtime,0);
    long time_usec=0;
    while ( testtime.tv_usec > updatetimeinmus) gettimeofday(&testtime,0);
    
    // main loop: give angle orders
    while(0==0){
        
        phase = iteration%numberofsteps;
        loop_number = floor(iteration/numberofsteps);
        
        if (loop_number%2==0) grid.high_duty = false; //to know when the pattern restarts to be read
        else grid.high_duty = true;
        
        // just read the values previously computed
        for (int i=0; i<numberOfServos;i++){
            positions[i]= stored_positions[i][phase];
            anglesteps[i]= stored_steps[i][phase];
        }
        
        //setposition of each servo:
        time_usec += updatetimeinmus;
        gettimeofday(&testtime,0);
        if(time_usec>1000000)time_usec-=1000000;
        if(testtime.tv_usec > time_usec) {
            cout << "---------------------------------Problem!!!------------------------------" << endl;
            cout << "difference: " << (time_usec - testtime.tv_usec) << " mu sec, testtime: " << testtime.tv_usec <<  " - time_usec: " << time_usec <<  endl;
        }
        while (testtime.tv_usec <= time_usec){
            gettimeofday(&testtime,0);
            if(time_usec==1000000 && testtime.tv_usec<updatetimeinmus ){
                break;
            }
        }
        setanglestoallservosII(positions,anglesteps,constant, target_rms); // for motion
        cout << iteration << "  " << phase << "  " << loop_number << "  " << grid.high_duty << "  \n";
        iteration++;
    }
    anglefile.close();
    return 0;
}
