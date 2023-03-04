#include "planets.h"

int main() {

    //declare values we will obtain from file
    double hf;
    double Gf;
    double Tf;
    int no_linesf = 0;

    //this part reads the file
    fstream vMyfile("parameters.txt", ios::in);

    //check that file is good
    if (vMyfile.good()) {
        
        bool line1 { true }; //this variable serves to indicate we are at the first line of the parameters file
        string line; //this variable is a string, containing the parameters file line.
        
        while (!vMyfile.eof()) { //read line by line the file to check for errors, and count the number of lines in the process
            
            ++no_linesf;
            
            try {

                getline(vMyfile, line, '\n');

                if (line.empty()) {
                    throw runtime_error("Parameters text file has empty lines, please delete these for code to run properly."); //possible error noticed when practicing
                }

                //split line into words
                stringstream ss(line);
                string datapoint;
              
                int count = 0;

                //count number of words
                while (ss >> datapoint) {
                    count++;
                } 
                
                //throw possible errors
                if (line1 && count!=3) {
                    throw runtime_error("The first line in the parameters file doesn't have exactly 3 prameters.");
                }

                if (!line1 && count != 5) {
                    throw runtime_error("At least one planet doesn't have exactly 5 parameters.");
                }
            }

            catch (const runtime_error& e) {
                cout << e.what() << endl;
                return 1;
            }

            line1 = false; //after the first iteration of the loop, we only have planet lines, so this can be set to false.
        }

        //after counting the number of lines go back to start to extract the data.
        vMyfile.seekg(0, ios::beg);
        getline(vMyfile, line);
        stringstream(line) >> Gf >> Tf >> hf;   

        //ensure variables remain constant through code
        const double h = hf;
        const double G = Gf;
        const double T = Tf;
        const int no_lines = no_linesf;

        //generate vector to hold all of our planet objects
        vector<Planet> planets(no_lines - 1);

        //assign values to each respective object from the file
        for (int i = 0; i < no_lines - 1; i++) {
            getline(vMyfile, line);
            double xtemp;
            double ytemp;
            double xdottemp;
            double ydottemp;
            double mtemp;
            stringstream(line) >> xtemp >> ytemp >> xdottemp >> ydottemp >> mtemp;    
            planets[i].SetX(xtemp); 
            planets[i].SetY(ytemp); 
            planets[i].SetXdot(xdottemp); 
            planets[i].SetYdot(ydottemp); 
            planets[i].SetM(mtemp);
        }

        //Close file, after finished operating on it
        vMyfile.close();

        //allow the output file to take in output
        ofstream vMyoutput("output.txt", ios::out);

        //Perform rk4 over the whole Time period, with the appropiate timestep
        for (double i = 0.0; i <= T+h; i += h) {
            //write to output file
            for (int z = 0; z < no_lines - 1; ++z) {
                vMyoutput << z + 1 << " " << i << " " << planets[z].GetX() << " " << planets[z].GetY() << " " << planets[z].GetXdot() << " " << planets[z].GetYdot() << endl;
            }
            //index 0 doesn't affect anything, just to call member function
            planets =  planets[0].rk4(planets, no_lines, G, h);
        }
        vMyoutput.close();
    }
    else {
        cout << "parameters.txt file not found. " << endl;
    }
} 