//
//  old methods.cpp
//  
//
//  Created by pilot_user on 9/17/12.
//
//



void algo::area(float actpos[], float actstep[]){
	//calculates the projected area and corrects the position of each servo to keep the projected area constant
	//float delta;
	float area_ratio;
	
	//create an array with servonumbers and shuffle it
	int random[numberOfServos];
	srand(time(0));
	for (int i=0; i<numberOfServos; i++) {
		random[i] = i;
	}
	for (int n=0; n < (numberOfServos - 1); n++) {
		int r = n + (rand() % (numberOfServos - n));
		int temp = random[n]; random[n] = random[r]; random[r] = temp;
	}
	
	
	projected_area = 0;
	if(angle_constant_area!=0){
	    for(int j=0;j<numberOfServos;j++){
			projected_area += 100 * fabs(sin(actpos[j] * M_PI / 180));
	    }
		//cout << "old projected area :" << projected_area << endl;
		//difference for each servo that needs to be corrected
		//cout << "constant area: " << angle_constant_area << endl;
		projected_area = (angle_constant_area - projected_area) / numberOfServos;
		//cout << "area per servo that needs to be corrected: " << projected_area << endl;
	    for(int l=0; l<numberOfServos; l++){
			area_ratio = (100 * fabs(sin(actpos[random[l]] * M_PI / 180)) + projected_area) / 100;
			if (area_ratio <= 1) {
				if (area_ratio > 0) {
					epsilon = (asin(area_ratio) * 180 / M_PI) - fabs(actpos[random[l]]);
				}
				else {
					epsilon = - fabs(actpos[random[l]]);
					area_ratio = area_ratio * 100;
					projected_area += area_ratio / (numberOfServos - l);
				}
		    }
			else {
				epsilon = 90 - fabs(actpos[random[l]]);
				area_ratio = (area_ratio * 100) - 100;
				projected_area += area_ratio / (numberOfServos - l);	
			}
			//cout << l << "  before:" << "actpos: " << actpos[l] << "/t" << "actstep " << actstep[l] << endl;
			//cout << "epsilon: " << epsilon << endl;
			overwrite(actpos, actstep, random[l]);
			//projected_area += delta / (numberOfServos - l);
		}
		projected_area = 0;
		for(int j=0;j<numberOfServos;j++){
			projected_area += 100 * fabs(sin(actpos[j] * M_PI / 180));
			//cout << j << "  " << "currently projected area: " << projected_area << endl;
			if (actstep[j] > 42.8 || actstep[j] < 0) {			
				cout << j << "  after:" << "actpos: " << actpos[j] << "/t" << "actstep " << actstep[j] << endl;
			}
		}
		cout << "new projected area: " << projected_area << endl;
	}
}


float algo::overwrite (float actpos[], float actstep[], int servonumber) {
	//corrects the arrays positions, steps, actpos and actstep
	//float check[positions[servonumber].size()];
	//int count =0;
	//float standard;
	//standard = 100 * fabs(sin((actpos[servonumber] + epsilon) * M_PI / 180));
	//cout << "standard alt: " << standard << endl;
	//while (true) {
	//for (int pos = actualpositioninvector[servonumber]; pos <= (int)(positions[servonumber].size() - 1); pos++) {
	//	check[pos] = positions[servonumber].at(pos);
	//}
	for (int vecpos = actualpositioninvector[servonumber] - 1; vecpos <= (int)(positions[servonumber].size() - 1); vecpos++) {
		if (positions[servonumber].at(vecpos) == -100) {
			break;
		}
		if (positions[servonumber].at(vecpos) < 0) {
			//if because of an calculated epsilon an angle is bigger than 90 or -90 degrees, this overstepping will move the
			//paddle in the other direction
			if ((positions[servonumber].at(vecpos) - epsilon)<= -90 ) {
				positions[servonumber].at(vecpos) = -180 + (epsilon - positions[servonumber].at(vecpos));
			}
			else positions[servonumber].at(vecpos) -= epsilon;
		}
		if (positions[servonumber].at(vecpos) >= 0) {
			if ((positions[servonumber].at(vecpos) + epsilon)>= 90) {
				positions[servonumber].at(vecpos) = 180 - (epsilon + positions[servonumber].at(vecpos));
			}
			else positions[servonumber].at(vecpos) += epsilon;
		}
	}	
	/*if (count == 1) {
	 for (int pos = 0; pos <= (int)(positions[servonumber].size() - 1); pos++) {
	 positions[servonumber].at(pos) = check[pos];
	 //formally here actpos[] = ...
	 }
	 standard = standard - (100 * fabs(sin(actpos[servonumber] * M_PI / 180)));
	 cout << "standard neu:" << standard << endl;
	 return standard;
	 break;
	 }*/
	
	//correction of the next anglestep
	if (actualpositioninvector[servonumber] - 1 == 0) {
		steps[servonumber].at(actualpositioninvector[servonumber] - 1) += epsilon;
	}
	else {
		if (positions[servonumber].at(actualpositioninvector[servonumber] - 1) > 0) {
			if (positions[servonumber].at(actualpositioninvector[servonumber] - 1) - positions[servonumber].at(actualpositioninvector[servonumber] -2) - (steps[servonumber].at(actualpositioninvector[servonumber] - 1) + epsilon) < 0.01 || positions[servonumber].at(actualpositioninvector[servonumber] -2) < 0) {
				//float difference = positions[servonumber].at(actualpositioninvector[servonumber] - 1) - positions[servonumber].at(actualpositioninvector[servonumber] -2) - (steps[servonumber].at(actualpositioninvector[servonumber] - 1) + epsilon);
				//cout << "actual position :" << positions[servonumber].at(actualpositioninvector[servonumber] - 1) << endl;
				//cout << "last position :" << positions[servonumber].at(actualpositioninvector[servonumber] - 2) << endl;
				//cout << "actual step :" << steps[servonumber].at(actualpositioninvector[servonumber] - 1) << endl;
				//cout << "epsilon :" << epsilon << endl;
				//cout << "difference :" << difference << endl;
				//steps[servonumber].at(actualpositioninvector[servonumber] -1) += epsilon;
			}
			else steps[servonumber].at(actualpositioninvector[servonumber] -1) -= epsilon;
		}
		if (positions[servonumber].at(actualpositioninvector[servonumber] - 1) < 0) {
			if(fabs(positions[servonumber].at(actualpositioninvector[servonumber] - 1) - positions[servonumber].at(actualpositioninvector[servonumber] -2)) - (steps[servonumber].at(actualpositioninvector[servonumber] -1) + epsilon) < 0.01 || positions[servonumber].at(actualpositioninvector[servonumber] -2) > 0) {
				steps[servonumber].at(actualpositioninvector[servonumber] -1) += epsilon;
			}
			else steps[servonumber].at(actualpositioninvector[servonumber] -1) -= epsilon;
		}
	}
	
	/*if (steps[servonumber].at(actualpositioninvector[servonumber] - 1) > max_angleperstep) {
	 if (epsilon > 0) {
	 epsilon -= (steps[servonumber].at(actualpositioninvector[servonumber] - 1) - max_angleperstep);
	 }
	 if (epsilon < 0) {
	 epsilon += (steps[servonumber].at(actualpositioninvector[servonumber] - 1) - max_angleperstep);
	 }
	 steps[servonumber].at(actualpositioninvector[servonumber] - 1) = max_angleperstep;
	 }
	 
	 if (steps[servonumber].at(actualpositioninvector[servonumber] - 1) < 0) {
	 if (epsilon > 0) {
	 epsilon -= fabs(steps[servonumber].at(actualpositioninvector[servonumber] - 1));
	 }
	 if (epsilon < 0) {
	 epsilon += fabs(steps[servonumber].at(actualpositioninvector[servonumber] - 1));
	 }
	 steps[servonumber].at(actualpositioninvector[servonumber] - 1) = 0;
	 }*/
	actstep[servonumber] = steps[servonumber].at(actualpositioninvector[servonumber] - 1);
	actpos[servonumber] = positions[servonumber].at(actualpositioninvector[servonumber] - 1);
	
	//count++;
	//}
}