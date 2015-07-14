/*  Program to print letters (and maybe .bmp if I have time to implement it) on the active grid.
This should not have a big scientific interest, but it was a nice thing to do.
 
 Called using ./Printing in terminal.
 
 For now there are just latin letters and a few special caracters. But do not hesitate to create the ones you need! Or the numbers!
 Fore each character, you just have to define the dark "pixels" (ie the closed wings).
 
 The way special characters are handled now is not satisfactory. If you have time to change that, it will be great! (I don't have time to plunge into the details of wilde characters handling in C++)
 
 Please enjoy! Viel Spaß! Amusez-vous bien ! :)
 
 Florent Lachaussée June 2013
 florent.lachaussee@ens.fr
*/


#include <iostream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <string>
#include <wchar.h>

using namespace std;

#include "lib/algo.h"

void write_letter (char letter);
int main (int argc , char * const argv[]);
void wait (float seconds);

double angles[14][12]; // to store the angles (used to control "pixel" color)
double inverted_angles[14][12];
algo alg; //just create an algo object


int main (int argc , char * const argv[]){

    char message[256]; // to store user's message
    int l =0; // to index the letters in message

    cout << endl << "Welcome to " << argv[0] << endl;
    cout << "This is the program to print things on the active grid  " << endl;
    cout << "Please type your message (use only capitals letters and '_' for spacing). \nSpecial caracters availables are: É (type &),Ü (type #) \n";
    cin >> message;

    cout << "You want to write: ";
    cout << message << endl;
    
    wait(10);
    
    char currentchar = message[0];
    
    while (currentchar != '\0')
    {   currentchar = message[l];
        cout << currentchar <<" \n";
        write_letter(currentchar);
        l++;
    }
    return 0;
}

void wait ( float seconds )
{
	clock_t endwait;
	endwait = clock () + seconds * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}

void write_letter(char letter){
    // set all paddles to 0 degrees
    for(int i=0; i<14; i++){
        for(int j=0; j<12; j++){
            angles[i][j]=0;
        }
    }
    
    
    // for each letter, set the interesting paddles to 90 degrees
    if(letter == 'A'){
        for(int i=5;i<10;i++){ // horizontal bar of the A
            angles[i][5]=90;
        }
        angles[4][3]=angles[4][4]=90;
        angles[5][6]=90;
        angles[6][7]=angles[6][8]=90;
        angles[7][8]=angles[7][9]=90;
        angles[8][7]=angles[8][8]=90;
        angles[9][6]=90;
        angles[10][3]=angles[10][4]=90;
    }
    
    else if(letter == 'B'){
        for(int j=3;j<10;j++){
            angles[6][j]=90;
        }
        angles[7][3]=angles[7][6]=angles[7][9]=90;
        angles[8][3]=angles[8][6]=angles[8][7]=angles[8][8]=90;
        angles[9][4]=angles[9][5]=90;
    }
    
    else if(letter == 'C'){
        for(int j=5;j<8;j++){
            angles[5][j]=90;
        }
        for(int i=7;i<10;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        angles[6][4]=angles[6][8]=90;

    }
    
    else if(letter == 'D'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        angles[6][3]=angles[6][9]=90;
        angles[7][3]=angles[7][9]=90;
        angles[8][4]=angles[8][8]=90;
        angles[9][5]=angles[9][6]=angles[9][7]=90;

    }
    
    else if(letter == 'E'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        for(int i=6;i<10;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        angles[6][6]=angles[7][6]=90;

    }
    
    else if(letter == 'F'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        for(int i=6;i<10;i++){
            angles[i][9]=90;
        }
        angles[6][6]=angles[7][6]=90;

    }
    
    else if(letter == 'G'){
        for(int j=5;j<8;j++){
            angles[5][j]=90;
        }
        for(int i=7;i<10;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        angles[6][4]=angles[6][8]=90;
        angles[8][5]=angles[9][4]=angles[9][5]=90;

    }
    
    else if(letter == 'H'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
            angles[9][j]=90;
        }
        for(int i=6;i<9;i++){
            angles[i][6]=90; // horizontal bar
        }

    }
    
    else if(letter == 'I'){
        for(int j=3;j<10;j++){
            angles[7][j]=90;
        }
        angles[6][3]=angles[6][9]=90;
        angles[8][3]=angles[8][9]=90;

    }
    
    else if(letter == 'J'){
        for(int j=3;j<10;j++){
            angles[7][j]=90;
        }
        for (int i=5; i<10;i++){
            angles[i][9]=90;
        }
        angles[4][3]=angles[4][4]=90;
        angles[5][2]=angles[6][2]=90;
    }

    else if(letter == 'K'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        for (int i=6; i<10;i++){
            angles[i][i]=90;
            angles[i][12-i]=90;
        }
    }
    
    else if(letter == 'L'){
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        for(int i=6;i<10;i++){
            angles[i][3]=90;
        }

    }
    
    else if(letter == 'M'){
        for(int j=3;j<10;j++){
            angles[4][j]=90;
            angles[10][j]=90;
        }
        angles[5][8]=angles[6][7]=90;
        angles[7][6]=90;
        angles[8][7]=angles[9][8]=90;

    }
    
    else if(letter == 'N'){
        for(int j=3;j<10;j++){
            angles[4][j]=90;
            angles[10][j]=90;
        }
        for(int i=5;i<10;i++){
            angles[i][13-i]=90; // diagonal
        }

    }
    
    else if(letter == 'O'){
        for(int j=5;j<8;j++){
            angles[4][j]=90;
            angles[10][j]=90;
        }
        for(int i=6;i<9;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        angles[5][4]=angles[5][8]=90;
        angles[9][4]=angles[9][8]=90;

    }
    
    else if(letter == 'P'){
        for(int j=3;j<10;j++){
            angles[6][j]=90;
        }
        angles[7][6]=angles[7][9]=90;
        angles[8][6]=angles[8][9]=90;
        angles[9][7]=angles[9][8]=90;

    }
    
    else if(letter == 'Q'){
        for(int j=5;j<8;j++){
            angles[4][j]=90;
            angles[10][j]=90;
        }
        for(int i=6;i<9;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        for(int i=8;i<11;i++){
            angles[i][13-i]=90;
        }
        angles[5][4]=angles[5][8]=90;
        angles[9][4]=angles[9][8]=90;
    }
    
    else if(letter == 'R'){
        for(int j=3;j<10;j++){
            angles[6][j]=90;
        }
        for (int i=7; i<10;i++){
            angles[i][12-i]=90;
        }
        angles[7][6]=angles[7][9]=90;
        angles[8][6]=angles[8][9]=90;
        angles[9][7]=angles[9][8]=90;

    }
    
    else if(letter == 'S'){
        for(int i=6;i<9;i++){
            angles[i][3]=90;
            angles[i][6]=90;
            angles[i][9]=90;
        }
        angles[5][3]=angles[5][7]=angles[5][8]=90;
        angles[9][4]=angles[9][5]=angles[9][9]=90;

    }
    
    else if(letter == 'T'){
        for(int j=3;j<10;j++){
            angles[7][j]=90;}
        for(int i=5;i<10;i++){
            angles[i][9]=90;
        }

    }
    
    else if(letter == 'U'){
        for(int j=4;j<10;j++){
            angles[5][j]=90;
            angles[9][j]=90;
        }
        for(int i=6;i<9;i++){
            angles[i][3]=90; // bottom
        }

    }
    
    else if(letter == 'V'){
        for(int i=4;i<7;i++){
            angles[i][16-2*i]=angles[i][17-2*i]=90;
        }
        for(int i=8;i<11;i++){
            angles[i][2*i-12]=angles[i][2*i-11]=90;
        }
        angles[7][3]=angles[7][4]=90;

    }
    
    else if(letter == 'W'){
        for(int i=2;i<5;i++){
            angles[i][13-2*i]=angles[i][12-2*i]=90;
        }
        for(int i=10;i<13;i++){
            angles[i][2*i-16]=angles[i][2*i-15]=90;
        }
        angles[5][3]=angles[5][4]=90;
        angles[6][4]=angles[6][5]=90;
        angles[7][6]=angles[7][7]=90;
        angles[8][4]=angles[8][5]=90;
        angles[9][3]=angles[9][4]=90;

    }
    
    else if(letter == 'X'){
        for(int i=4;i<11;i++){
            angles[i][13-i]=90;
            angles[i][i-1]=90;
        }
    }
    
    else if(letter == 'Y'){
        for(int j=2;j<7;j++){
            angles[7][j]=90;
        }
        angles[4][9]=angles[5][8]=angles[6][7]=90;
        angles[8][7]=angles[9][8]=angles[10][9]=90;
    }
    
    else if(letter == 'Z'){
        for(int i=4;i<11;i++){
            angles[i][3]=90;
            angles[i][9]=90;
            angles[i][i-1]=90;
        }
        
    }
    
    else if(letter == '&'){// E accent aigu
        for(int j=3;j<10;j++){
            angles[5][j]=90;
        }
        for(int i=6;i<10;i++){
            angles[i][3]=90;
            angles[i][9]=90;
        }
        angles[6][6]=angles[7][6]=90;
        
        angles[7][10]=angles[8][11]=90; // accent aigu
    }

    else if(letter == '#'){ // U Umlaut
        for(int j=4;j<10;j++){
            angles[5][j]=90;
            angles[9][j]=90;
        }
        for(int i=6;i<9;i++){
            angles[i][3]=90; // bottom
        }
        angles[6][10]=angles[6][11]=90; // Umlaut
        angles[8][10]=angles[8][11]=90; // Umlaut
    }
    
    else if(letter == '_'){ // open the grid
    }
    
    else {cout << "NON!!" << endl; // Non possible character
} // The grid will open (as for blank)
    
    // send the orders to the servos
    alg.grid.setanglesII(angles);
    
    // inversion to have better contrast
    /*for(int i=0; i<14; i++){
        for(int j=0; j<12; j++){
            if (angles[i][j]==0)
            {inverted_angles[i][j]=90;}
            else if (angles[i][j]==90)
            {inverted_angles[i][j]=0;}
        }
    }
    alg.grid.setanglesII(inverted_angles);*/
    
    wait(1); // wait in order to be able to see the letter a long time enough on the grid
}