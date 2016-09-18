//Program for Project
//Vashist Kutti Suresh
//Email: vsuresh@clemson.edu
//Date : 12/3/2015

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <OpenImageIO/imageio.h>

OIIO_NAMESPACE_USING

#include "Matrix.h"
#include "Vector.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std;
static int WIDTH,WIDTH2,FWIDTH,SWIDTH,PWIDTH;
static int HEIGHT,HEIGHT2,FHEIGHT,SHEIGHT,PHEIGHT;
static int CHANNELS,FCHANNELS,SCHANNELS,PCHANNELS;
static char INPUT_FILE[512];
static char Write_file[512];
float leftpix,rightpix,toppix,bottompix;
float xc,yc,strength,centerx,centery,minDim;
float scalex,scaley;
static float **lw;
static float **ld;
int orientationchoice;
int channelchoice;
static int gammaflag = 0;
static int redflag = 0;
static int greenflag = 0;
static int blueflag = 0;
static int hueflag = 0;
static int saturationflag = 0;
static int valueflag = 0;
static int compareflag = 1;
int choice;

class rgba_pixel {
   public:
      float r;
      float g;
      float b;
      float a;

};

static rgba_pixel ** undo_buffer;
static rgba_pixel ** original_buffer;
static rgba_pixel ** parisflag_buffer;
static rgba_pixel ** flare_buffer;
static rgba_pixel ** frame_buffer;
static rgba_pixel ** pix_buffer2;
static rgba_pixel ** output_buffer;

static ImageSpec the_spec;
static vector<float> pixels;
static vector<float> the_oiio_pixels;
static vector<float> the_oiio_pixels1;
static vector<float> the_oiio_pixels2;
static vector<float> the_oiio_pixels3;

static float r = 1.0;
static float g = 1.0;
static float b = 1.0;
static float h = 0.0;
static float s = 0.0;
static float v = 0.0;

Matrix3x3 M(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
Matrix3x3 I(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

//Function to copy the output buffer to the original buffer , for using it while comparing  
void copyoriginal(){
	
	original_buffer = new rgba_pixel*[HEIGHT];
	original_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      original_buffer[i] = original_buffer[i-1] + WIDTH;
	}
	
	for(int row = 0; row < HEIGHT; row++){
		for(int col = 0; col < WIDTH; col++) {
			original_buffer[row][col].r = output_buffer[row][col].r;
			original_buffer[row][col].g = output_buffer[row][col].g;
			original_buffer[row][col].b = output_buffer[row][col].b;
			original_buffer[row][col].a = 1.0;
			}
		}
}

//Funtion to backup the output_buffer into the undo_buffer
void backup(){
	
	undo_buffer = new rgba_pixel*[HEIGHT];
	undo_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      undo_buffer[i] = undo_buffer[i-1] + WIDTH;
	}
	
	for(int row = 0; row < HEIGHT; row++){
		for(int col = 0; col < WIDTH; col++) {
			undo_buffer[row][col].r = output_buffer[row][col].r;
			undo_buffer[row][col].g = output_buffer[row][col].g;
			undo_buffer[row][col].b = output_buffer[row][col].b;
			undo_buffer[row][col].a = 1.0;
			}
		}
}
//Function to take backup of the last action performed output buffer
void restorebackup(){
	
	for(int row = 0; row < HEIGHT; row++){
		for(int col = 0; col < WIDTH; col++) {
			output_buffer[row][col].r = undo_buffer[row][col].r;
			output_buffer[row][col].g = undo_buffer[row][col].g;
			output_buffer[row][col].b = undo_buffer[row][col].b;
			output_buffer[row][col].a = 1.0;
			}
		}
}

//Function to populate opengl output buffer
void populate_output_buffer() {                                
   
   for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) {
         output_buffer[HEIGHT-row-1][col].r = the_oiio_pixels[(row*WIDTH+col)*CHANNELS + 0];
         output_buffer[HEIGHT-row-1][col].g = the_oiio_pixels[(row*WIDTH+col)*CHANNELS + 1];
         output_buffer[HEIGHT-row-1][col].b = the_oiio_pixels[(row*WIDTH+col)*CHANNELS + 2];
          if (CHANNELS == 4) {
            output_buffer[HEIGHT-row-1][col].a = the_oiio_pixels[(row*WIDTH+col)*CHANNELS + 3];
         } else {
            output_buffer[HEIGHT-row-1][col].a = 1.0;                 
         }
      }         
      }
   CHANNELS = 4;
   copyoriginal();
}

//Function to populate the frame buffer after loading the frame image
void populate_frame_buffer() { 
                            
   for (int row=0; row<FHEIGHT; row++) {
      for (int col=0; col<FWIDTH; col++) {

         frame_buffer[FHEIGHT-row-1][col].r = the_oiio_pixels2[(row*FWIDTH+col)*FCHANNELS + 0];
         frame_buffer[FHEIGHT-row-1][col].g = the_oiio_pixels2[(row*FWIDTH+col)*FCHANNELS + 1];
         frame_buffer[FHEIGHT-row-1][col].b = the_oiio_pixels2[(row*FWIDTH+col)*FCHANNELS + 2];
          if (FCHANNELS == 4) {
            frame_buffer[FHEIGHT-row-1][col].a = the_oiio_pixels2[(row*FWIDTH+col)*FCHANNELS + 3];
         } else {
           
            frame_buffer[FHEIGHT-row-1][col].a = 1.0;                  
         }
      }        
      }
   
}

//Function to populate the flare buffer after loading the flare image
void populate_flare_buffer() { 
                              
   for (int row=0; row<SHEIGHT; row++) {
      for (int col=0; col<SWIDTH; col++) {

         flare_buffer[SHEIGHT-row-1][col].r = the_oiio_pixels1[(row*SWIDTH+col)*SCHANNELS + 0];
         flare_buffer[SHEIGHT-row-1][col].g = the_oiio_pixels1[(row*SWIDTH+col)*SCHANNELS + 1];
         flare_buffer[SHEIGHT-row-1][col].b = the_oiio_pixels1[(row*SWIDTH+col)*SCHANNELS + 2];
          if (SCHANNELS == 4) {
            flare_buffer[SHEIGHT-row-1][col].a = the_oiio_pixels1[(row*SWIDTH+col)*SCHANNELS + 3];
         } else {
           
            flare_buffer[SHEIGHT-row-1][col].a = 1.0;                 
         }
       }
     }
   
}

//Funtion to populate the paris flag buffer after loading the paris flag image
void populate_parisflag_buffer() { 
                 
	parisflag_buffer = new rgba_pixel*[HEIGHT];
	parisflag_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      parisflag_buffer[i] = parisflag_buffer[i-1] + WIDTH;
	}  
	 
   for (int row=0; row<PHEIGHT; row++) {
      for (int col=0; col<PWIDTH; col++) {

         parisflag_buffer[PHEIGHT-row-1][col].r = the_oiio_pixels3[(row*PWIDTH+col)*PCHANNELS + 0];
         parisflag_buffer[PHEIGHT-row-1][col].g = the_oiio_pixels3[(row*PWIDTH+col)*PCHANNELS + 1];
         parisflag_buffer[PHEIGHT-row-1][col].b = the_oiio_pixels3[(row*PWIDTH+col)*PCHANNELS + 2];
          if (PCHANNELS == 4) {
            parisflag_buffer[PHEIGHT-row-1][col].a = the_oiio_pixels3[(row*PWIDTH+col)*PCHANNELS + 3];
         } else {
           
            parisflag_buffer[PHEIGHT-row-1][col].a = 1.0;                  
         }
      }       
    }
   
}


//Function to load the image which needs to be edited
void loadImage(char* input_file){                       
	
	ImageInput *in = ImageInput::open(input_file);
	if (!in) {
		cerr << "could not open file" << endl;
		exit(-1);
	}

	the_spec = in->spec();
	WIDTH = the_spec.width;
	HEIGHT = the_spec.height;
	CHANNELS = the_spec.nchannels;
	the_oiio_pixels.resize(WIDTH*HEIGHT*CHANNELS);
	in->read_image (TypeDesc::FLOAT, &the_oiio_pixels[0]);

	in->close ();
	delete in;
}

//Function to load the frame image which is applied on the editing image
void loadFrame(char* input_file){                       
	
	ImageInput *in = ImageInput::open(input_file);
	if (!in) {
		cerr << "could not open file" << endl;
		exit(-1);
	}

	the_spec = in->spec();
	FWIDTH = the_spec.width;
	FHEIGHT = the_spec.height;
	FCHANNELS = the_spec.nchannels;	
	the_oiio_pixels2.resize(FWIDTH*FHEIGHT*FCHANNELS);
	in->read_image (TypeDesc::FLOAT, &the_oiio_pixels2[0]);
	populate_frame_buffer();
	in->close ();
	delete in;
}

//Function to load the flare image which is to be applied on the editing image
void loadFlare(char* input_file){                  
	
	ImageInput *in = ImageInput::open(input_file);
	if (!in) {
		cerr << "could not open file" << endl;
		exit(-1);
	}

	the_spec = in->spec();
	SWIDTH = the_spec.width;
	SHEIGHT = the_spec.height;
	SCHANNELS = the_spec.nchannels;
	the_oiio_pixels1.resize(SWIDTH*SHEIGHT*SCHANNELS);
	in->read_image (TypeDesc::FLOAT, &the_oiio_pixels1[0]);
	populate_flare_buffer();
	in->close ();
	delete in;
}

//Function to load the frame image which is applied on the editing image
void loadParis(char* input_file){                    
	ImageInput *in = ImageInput::open(input_file);
	if (!in) {
		cerr << "could not open file" << endl;
		exit(-1);
	}

	the_spec = in->spec();
	PWIDTH = the_spec.width;
	PHEIGHT = the_spec.height;
	PCHANNELS = the_spec.nchannels;	
	the_oiio_pixels3.resize(PWIDTH*PHEIGHT*PCHANNELS);
	in->read_image (TypeDesc::FLOAT, &the_oiio_pixels3[0]);
	populate_parisflag_buffer();
	in->close ();
	delete in;
}


//Function to write the image
int writeimg(const char * Write_file){
	
	if(Write_file == NULL) {
    cout<<"Write filename not provided."<<endl;
    exit(-1);
  }  
  
  ImageOutput * out = ImageOutput::create(Write_file);
  if (!out) {
    std::cerr << "Write file could not created " << Write_file << geterror();
  }

  ImageSpec spec(WIDTH, HEIGHT, CHANNELS, TypeDesc::FLOAT);
  
  out->open(Write_file , spec);
  pixels.resize(WIDTH*HEIGHT*CHANNELS);
  int i = 0;
  for(int row = HEIGHT-1; row >= 0; --row) {
    for(int col = 0; col < WIDTH; ++ col) {
      pixels[i++] = output_buffer[row][col].r;
      pixels[i++] = output_buffer[row][col].g;
      pixels[i++] = output_buffer[row][col].b;
      if(CHANNELS == 4)
        pixels[i++] = output_buffer[row][col].a;
    }
  }

  if(out->write_image(TypeDesc::FLOAT, &pixels[0]))
    cout<<"Image saved"<<endl;
  else
    cout<<"Write Error "<<geterror()<<endl;
  out->close();
  delete out;

}

//Function to draw and display the image on the window using glDrawPixels()
void drawImage(){     	                                      
	glClear(GL_COLOR_BUFFER_BIT); 
	glViewport(0,0,WIDTH,HEIGHT);
	glRasterPos2i(0,0); 
	glDrawPixels(WIDTH,HEIGHT,GL_RGBA,GL_FLOAT,output_buffer[0]);
	glFlush();
}

//Funtion to draw the original image while comparing operation is called 
void drawOriginal(){     	                                      
		
	glClear(GL_COLOR_BUFFER_BIT); 
	glViewport(0,0,WIDTH,HEIGHT);
	glRasterPos2i(0,0); 
	glDrawPixels(WIDTH,HEIGHT,GL_RGBA,GL_FLOAT,original_buffer[0]);
	glFlush();
}

//Funtion to calculate the minimum among 4 values
int minimum(int v1,int v2,int v3,int v4){
    
    int min1,min2,minm;
  
	min1 = min(v1,v2);
	min2 = min(v3,v4);
	minm = min(min1,min2);
			
	return minm;
}
  
//Funtion to calculate the maximum among 4 values
int maximum(int v1,int v2,int v3,int v4){
    
    int max1,max2,maxm;
		
	max1 = max(v1,v2);
	max2 = max(v3,v4);
	maxm = max(max1,max2);
			
	return maxm;
}

//Function to calculate the inverse of the matrix and 
//perform allocation of pixel using Inverse transform
void inversetransform(Matrix3x3 &M){
	
	Matrix3x3 I(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
//Compute the inverse of the matrix M	
	I = M.inv();
	
		Vector3d origin(leftpix, bottompix, 0);
	for (int y = 0; y < HEIGHT2; y++) {
		for (int x = 0; x < WIDTH2; x++) { 
			Vector3d pixel_out(x, y, 1); 
			pixel_out = pixel_out + origin;
			Vector3d pixel_in = I * pixel_out;
         
			float u = pixel_in[0] / pixel_in[2];
			float v = pixel_in[1] / pixel_in[2];
    
			int u1 = (int)round(u);
			int v1 = (int)round(v);

//Checking the boundary conditions of the image and giving 0.0 for outerbounds       
			if(u1<0 || v1<0 || u1>=WIDTH || v1>=HEIGHT)
			{
				pix_buffer2[y][x].r = 0.0;
				pix_buffer2[y][x].g = 0.0;
				pix_buffer2[y][x].b = 0.0;
				pix_buffer2[y][x].a = 1.0;
			}
			else
			{
				pix_buffer2[y][x].r = output_buffer[v1][u1].r;
				pix_buffer2[y][x].g = output_buffer[v1][u1].g;
				pix_buffer2[y][x].b = output_buffer[v1][u1].b;
				pix_buffer2[y][x].a = 1.0;
			}    
		}
	}

	HEIGHT=HEIGHT2;
	WIDTH=WIDTH2;

   output_buffer = new rgba_pixel*[HEIGHT];
   output_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
   for (int i=1; i<HEIGHT; i++) {
      output_buffer[i] = output_buffer[i-1] + WIDTH;
   } 
	for (int y = 0; y < HEIGHT; y++) {
			for (int x = 0; x < WIDTH; x++) { 
				output_buffer[y][x].r = pix_buffer2[y][x].r;
				output_buffer[y][x].g = pix_buffer2[y][x].g;
				output_buffer[y][x].b = pix_buffer2[y][x].b;
				output_buffer[y][x].a = pix_buffer2[y][x].a;
			}
		}
}


//Funtion to perform forward transform 
//get the left,right,top,bottom pixel , using which calculate new Height and Width
//and allocate new pixel buffer with the height and width 
void forwardtransform(Matrix3x3 &M){
	
	int h = HEIGHT;
	int w = WIDTH;
	Vector3d v1(0, 0, 1);
	Vector3d v2(0, h-1, 1);
	Vector3d v3(w-1, 0, 1);
	Vector3d v4(w-1, h-1, 1);

//Multiplying the Accumulated matrix with the vectors denoting the corners		
	v1 = M * v1;
	v1[0] = v1[0]/v1[2];
	v1[1] = v1[1]/v1[2];
	v1[2] = v1[2]/v1[2];
		 		 
	v2 = M * v2;
	v2[0] = v2[0]/v2[2];
	v2[1] = v2[1]/v2[2];
	v2[2] = v2[2]/v2[2];
		
	v3 = M * v3;
	v3[0] = v3[0]/v3[2];
	v3[1] = v3[1]/v3[2];
	v3[2] = v3[2]/v3[2];
		
	v4 = M * v4;
	v4[0] = v4[0]/v4[2];
	v4[1] = v4[1]/v4[2];
	v4[2] = v4[2]/v4[2];

//Determining the left,right,top and bottom pixel		
	leftpix = minimum(v1[0],v2[0],v3[0],v4[0]);
	rightpix = maximum(v1[0],v2[0],v3[0],v4[0]);
	bottompix = minimum(v1[1],v2[1],v3[1],v4[1]);
	toppix = maximum(v1[1],v2[1],v3[1],v4[1]);	  

//Calculating the new Width and Height					  
	HEIGHT2 = toppix-bottompix+1;
	WIDTH2 = rightpix-leftpix+1;
	
//Allocating and intializing the new pixelbuffer   
   pix_buffer2 = new rgba_pixel*[HEIGHT2];
   pix_buffer2[0] = new rgba_pixel[WIDTH2*HEIGHT2];
   for (int i=1; i<HEIGHT2; i++) {
      pix_buffer2[i] = pix_buffer2[i-1] + WIDTH2;
   } 
 
//Calling for the computation of the inverse transform   
	inversetransform(M);
}


//Function to generate the transformation matrix of rotation
// Multiply M by a rotation matrix of angle theta
void Rotate(Matrix3x3 &M, float theta){
  
   int row, col;
   Matrix3x3 R(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
   double rad, c, s;
   
   rad = PI * theta / 180.0;
   c = cos(rad);
   s = sin(rad);
   
   R[0][0] = c;
   R[0][1] = -s;
   R[1][0] = s;
   R[1][1] = c;

   Matrix3x3 Prod = R * M;

   for(row = 0; row < 3; row++) {
      for(col = 0; col < 3; col++) {
         M[row][col] = Prod[row][col];
      }
   }

//Calling forwardtransform function with the transformation matrix as parameter   
   forwardtransform(M);
}

//Function to generate the transformation matrix of Scaling
// Multiply M by a scaling matrix of scale x factor and scale y factor
void Scale(Matrix3x3 &M, float scalex, float scaley){
  
   int row, col;
   Matrix3x3 S(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
   float x,y;
   x = scalex;
   y = scaley;
   	
   S[0][0] = x;
   S[1][1] = y;

   Matrix3x3 Prod = S * M;

   for(row = 0; row < 3; row++) {
      for(col = 0; col < 3; col++) {
         M[row][col] = Prod[row][col];
      }
   }
   
//Calling forwardtransform function with the transformation matrix as parameter    
   forwardtransform(M);
}

//Function to composite the frame on the output buffer 
//which contains the last edited image data
void framecomposite(){

	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) {  
		   output_buffer[HEIGHT-row-1][col].r = (frame_buffer[HEIGHT-row-1][col].a) * (frame_buffer[HEIGHT-row-1][col].r) + (1 - (frame_buffer[HEIGHT-row-1][col].a)) * (output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (frame_buffer[HEIGHT-row-1][col].a) * (frame_buffer[HEIGHT-row-1][col].g) + (1 - (frame_buffer[HEIGHT-row-1][col].a)) * (output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (frame_buffer[HEIGHT-row-1][col].a) * (frame_buffer[HEIGHT-row-1][col].b) + (1 - (frame_buffer[HEIGHT-row-1][col].a)) * (output_buffer[HEIGHT-row-1][col].b);
		   
		   output_buffer[HEIGHT-row-1][col].a = 1.0;
		}
	}
}

//Function to perform the overlay effect on the image
void overlay(){
	
	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) {
		   output_buffer[HEIGHT-row-1][col].r = (0.25) * (output_buffer[HEIGHT-row-1][col].r) + (1) * (output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (0.25) * (output_buffer[HEIGHT-row-1][col].g) + (1) * (output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (0.25) * (output_buffer[HEIGHT-row-1][col].b) + (1) * (output_buffer[HEIGHT-row-1][col].b);
		   
		   output_buffer[HEIGHT-row-1][col].a = 1;
		}
	}
}

//Function to apply the warming and cooling filter to the image
void wceffects(){	
	
	cout<<" Warmth and cooling filter"<<endl;
	cout<<" Select any one (use number specified) : 1.warm  2.cooling "<<endl;
	int wcchoice;
	scanf("%d",&wcchoice);
	if(wcchoice==1)
	{
	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) { 
		   output_buffer[HEIGHT-row-1][col].r = (0.25) * ((1)*output_buffer[HEIGHT-row-1][col].r) + (1) * (output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (0.25) * ((0.6)*output_buffer[HEIGHT-row-1][col].g) + (1) * (output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (0.25) * (0*output_buffer[HEIGHT-row-1][col].b) + (1) * (output_buffer[HEIGHT-row-1][col].b);
		   
		   output_buffer[HEIGHT-row-1][col].a = 1;
			}
		}
	}
	else if(wcchoice==2)
	{
		for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) { 
		   output_buffer[HEIGHT-row-1][col].r = (0.25) * ((0.2)*output_buffer[HEIGHT-row-1][col].r) + (1) * (output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (0.25) * ((0.8)*output_buffer[HEIGHT-row-1][col].g) + (1) * (output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (0.25) * ((0.9)*output_buffer[HEIGHT-row-1][col].b) + (1) * (output_buffer[HEIGHT-row-1][col].b);
		   
		   output_buffer[HEIGHT-row-1][col].a = 1;
			}
		}
	}
	drawImage();
}

//Function to get the frame choice from user,call loadFrame and apply the frame
void setframes(){
	
	int framechoice;
	cout<<" Select the type of frame from the choice 1.Stick frame 2.Photo frame 3.Grunge frame 4.Film 5.Grunge2 frame "<<endl;
	scanf("%d", &framechoice);
	char* frame[5] = {"one.png","two.png","three.png","four.png","five.png"};	
	loadFrame(frame[framechoice-1]);		
	framecomposite();
    cout<<"Done settings frames"<<endl;
}

//Function to perform horizontal, vertical flip and call the rotate function
void Orientation(int type){	
	
	float theta; 
   
	if (type == 1)
 	{	
		cout<<WIDTH<<endl;
		cout<<"Horizontal flip performed"<<endl;
		for(int row = 0; row < HEIGHT; row++){
			for(int col = 0; col < WIDTH; col++) {
				output_buffer[row][col].r = undo_buffer[row][WIDTH-col-1].r;
				output_buffer[row][col].g = undo_buffer[row][WIDTH-col-1].g;
				output_buffer[row][col].b = undo_buffer[row][WIDTH-col-1].b;
				output_buffer[row][col].a = 1.0;
			}
		}
	}
	else if (type == 2)
	{ 
		drawImage();
		cout<<"Vertical flip performed"<<endl;
		for(int row = 0; row < HEIGHT; row++){
			for(int col = 0; col < WIDTH; col++) {
				output_buffer[row][col].r = undo_buffer[HEIGHT-row-1][col].r;
				output_buffer[row][col].g = undo_buffer[HEIGHT-row-1][col].g;
				output_buffer[row][col].b = undo_buffer[HEIGHT-row-1][col].b;
				output_buffer[row][col].a = 1.0;
				}
			}
	}
	else if(type == 3)
	{
		cout<<"Rotation has been selected ,Enter angle of rotation"<<endl;
		scanf("%f",&theta);
		if(theta==0)
			cout<<"error angle"<<endl;
		else	
			Rotate(M,theta);
	}
}

//Function to perform RGB to HSV conversion
void RGBtoHSV (float r, float g, float b, float &h, float &s, float &v) {

	float red, green, blue;
	float maxc, minc, delta;

// r, g, b to 0 - 1 scale
	red = r ; green = g ; blue = b ;  
	maxc = max(max(red, green), blue);
	minc = min(min(red, green), blue);
// value is maximum of r, g, b
	v = maxc;       
// saturation and hue 0 if value is 0
	if(maxc == 0){    
		s = 0;
		h = 0;
	} else {
		s = (maxc - minc) / maxc; 	

		delta = maxc - minc;
		if(delta == 0)           
			h = 0;
		else{
			if(red == maxc)       
				h = (green - blue) / delta;
			else if(green == maxc)
				h = 2.0 + (blue - red) / delta;
			else 
				h = 4.0 + (red - green) / delta;
			h = h * 60.0;
			if(h < 0)
				h = h + 360.0;
		}
	}
}

//Function to perform HSV to RGB conversion
void HSVtoRGB ( float h, float s, float v, float &r, float &g, float &b ){
	
	int i;
	float f, p, q, t, red, green, blue;

	if (s == 0) {
		red = green = blue = v;
	} else {
		h /= 60.0;
		i = (int) floor(h);
		f = h - i;
		p = v * (1-s);
		q = v * (1-s*f);
		t = v * (1 - s * (1 - f));

		switch (i) {
			case 0:
				red = v;
				green = t;
				blue = p;
				break;
			case 1:
				red = q;
				green = v;
				blue = p;
				break;
			case 2:
				red = p;
				green = v;
				blue = t;
				break;
			case 3:
				red = p;
				green = q;
				blue = v;
				break;
			case 4:
				red = t;
				green = p;
				blue = v;
				break;
			default:
				red = v;
				green = p;
				blue = q;
				break;
		}
	}

	r = red ;
	g = green;
	b = blue;	
}

//Function to adjust the Hue,Saturation and value based on the user selection
void HSVcompute(int hsv,int direction){   
	
	for (int y = 0; y < HEIGHT; y++) {
		for (int x = 0; x < WIDTH; x++) { 
			
			float r1 = output_buffer[y][x].r;
			float g1 = output_buffer[y][x].g;
			float b1 = output_buffer[y][x].b;
//Converting the pixels from RGB to HSV color space			
			RGBtoHSV(r1,g1,b1,h,s,v);
			
			if(direction==1){
				if(hsv==1){
					h=h+10;
					hueflag = 1;
					saturationflag = 0;
					valueflag = 0;
				}
				else if(hsv==2){
					s=s+0.01;
					hueflag = 0;
					saturationflag = 1;
					valueflag = 0;
				}
				else if(hsv==3){
					v=v+0.05;
					hueflag = 0;
					saturationflag = 0;
					valueflag = 1;
				}
			}
			else if(direction==0){
				if(hsv==1){
					h=h-10;
					hueflag = 1;
					saturationflag = 0;
					valueflag = 0;
				}
				else if(hsv==2){
					s=s-0.01;
					hueflag = 0;
					saturationflag = 1;
					valueflag = 0;
				}
				else if(hsv==3){
					v=v-0.05;
					hueflag = 0;
					saturationflag = 0;
					valueflag = 1;
				}
				}
				else{
					cout<<"no change"<<endl;
				}					
	
			float r3,g3,b3;

//Converting the pixels back to RGB space			
			HSVtoRGB(h,s,v,r3,g3,b3);
			output_buffer[y][x].r = r3;
			output_buffer[y][x].g = g3;
			output_buffer[y][x].b = b3;
			output_buffer[y][x].a = output_buffer[y][x].a;
		}
	}
	drawImage();
}

//Function to select which channel to modify in HSV mode
void HSV(){
	
	cout<<"select which channel you want to increase/decrease 1.Hue 2.Saturation 3.Value"<<endl;
	int hsvchoice;
	scanf("%d",&hsvchoice);
	if(hsvchoice==1)
	{
		hueflag = 1;
		saturationflag = 0;
		valueflag = 0;
		cout<<"HUE mode activated"<<endl;
		cout<<"Press UP & DOWN keys to increase or decrease"<<endl;
	}
	else if(hsvchoice==2)
	{
		hueflag = 0;
		saturationflag = 1;
		valueflag = 0;
		cout<<"SATURATION mode activated"<<endl;
		cout<<"Press UP & DOWN keys to increase or decrease"<<endl;
	}
	else if(hsvchoice==3)
	{
		hueflag = 0;
		saturationflag = 0;
		valueflag = 1;
		cout<<"VALUE mode activated"<<endl;
		cout<<"Press UP & DOWN keys to increase or decrease"<<endl;
	}
	else
	{
		cout<<"invalid hsv choice"<<endl;
	}
}

//Function to perform gamma correction on the image
void gammacorrection(int dir){
	
	float gama = 0.5;
	cout<<"Press up and down arrows to adjust gamma values"<<endl;
	gammaflag = 1;
	
	if(dir==1){
		gama = gama+0.01;
	}
	else if(dir==0){
		gama = -gama-0.01;
	}
	else{
	cout<<""<<endl;
	}
	cout<<gama<<endl;
	lw = new float*[HEIGHT];
    lw[0] = new float[WIDTH*HEIGHT];   
    for (int i=1; i<HEIGHT; i++) 
    {
      lw[i] = lw[i-1] + WIDTH;
    }
    ld = new float*[HEIGHT];
    ld[0] = new float[WIDTH*HEIGHT];
    for (int i=1; i<HEIGHT; i++) {
      ld[i] = ld[i-1] + WIDTH;
    }
   
	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) 
      {
		  
		   lw[row][col] = ( ( 20.0 * output_buffer[row][col].r ) + ( 40.0 * output_buffer[row][col].g ) + output_buffer[row][col].b)/61.0;			 
		   ld[row][col] = exp(gama * log10(lw[row][col]));
		   
		   if(lw[row][col] < 0 )
		   {
				lw[row][col] = FLT_MIN;
		   }
			
			output_buffer[row][col].r =  ld[row][col] * output_buffer[row][col].r;
			output_buffer[row][col].g =  ld[row][col] * output_buffer[row][col].g;
			output_buffer[row][col].b =  ld[row][col] * output_buffer[row][col].b;
			output_buffer[row][col].a =	 1.0;
      }
   }
}

//Function to perform the mirror effect on the image
void mirroreffect(){
	
	int mirrorchoice;
	cout<<"select the tyoe of mirror effect 1 or 2"<<endl;
	scanf("%d",&mirrorchoice);
	if (mirrorchoice == 1){
	cout<<"Mirror effect 1"<<endl;
	for(int row = 0; row < HEIGHT; row++){
		for(int col = 0; col < WIDTH/2; col++) {
			output_buffer[row][col].r = output_buffer[row][WIDTH-col-1].r;
			output_buffer[row][col].g = output_buffer[row][WIDTH-col-1].g;
			output_buffer[row][col].b = output_buffer[row][WIDTH-col-1].b;
			output_buffer[row][col].a = 1.0;
			}
		}
	}
	else if(mirrorchoice==2){
	cout<<"Mirror effect 2"<<endl;
	for(int row = 0; row < HEIGHT/2; row++){
		for(int col = 0; col < WIDTH; col++) {
			output_buffer[row][col].r = output_buffer[HEIGHT-row-1][col].r;
			output_buffer[row][col].g = output_buffer[HEIGHT-row-1][col].g;
			output_buffer[row][col].b = output_buffer[HEIGHT-row-1][col].b;
			output_buffer[row][col].a = 1.0;
			}
		}
	}
	else
	{
		cout<<"Error mirror effect choice"<<endl;
	}
}

//Function to add flare onto the image
void addflare(){	
	
	cout<<"Flare effect selected"<<endl;
	for (int row=0; row<HEIGHT; row++) {
		for (int col=0; col<WIDTH; col++) {
		   output_buffer[HEIGHT-row-1][col].r = (flare_buffer[HEIGHT-row-1][col].a) * ((1.25)*flare_buffer[HEIGHT-row-1][col].r) + (1 - (flare_buffer[HEIGHT-row-1][col].a)) * ((1)*output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (flare_buffer[HEIGHT-row-1][col].a) * (flare_buffer[HEIGHT-row-1][col].g) + (1 - (flare_buffer[HEIGHT-row-1][col].a)) * ((1)*output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (flare_buffer[HEIGHT-row-1][col].a) * (flare_buffer[HEIGHT-row-1][col].b) + (1 - (flare_buffer[HEIGHT-row-1][col].a)) * ((0.5)*output_buffer[HEIGHT-row-1][col].b);
		   output_buffer[HEIGHT-row-1][col].a = 1.0;
		}
	}	
}

//Function to perform the Deserted effect
void sunshine(){
	
	int flarechoice;
	cout<<" Select the type of Deserted effect from the choice 1,2 "<<endl;
	scanf("%d", &flarechoice);
	char* flare[2] = {"leftsun.png","rightsun.png"};	
	loadFlare(flare[flarechoice-1]);		
	addflare();
}

//Funciton to perform the Pray for paris effect
//The paris flag is composited in the image with half alpha vaule
void paris(){
	
	loadParis("france.png");
	cout<<"Pray for paris flag filter"<<endl;

	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) { 
		   output_buffer[HEIGHT-row-1][col].r = (0.5) * (parisflag_buffer[HEIGHT-row-1][col].r) + (0.5) * (output_buffer[HEIGHT-row-1][col].r);
		   output_buffer[HEIGHT-row-1][col].g = (0.5) * (parisflag_buffer[HEIGHT-row-1][col].g) + (0.5) * (output_buffer[HEIGHT-row-1][col].g);
		   output_buffer[HEIGHT-row-1][col].b = (0.5) * (parisflag_buffer[HEIGHT-row-1][col].b) + (0.5) * (output_buffer[HEIGHT-row-1][col].b);		   
		   output_buffer[HEIGHT-row-1][col].a = 1.0;
		}
	}	
}

//Function to compute a greyscale image ,when the user selects the Black and white effect
void greyscaleeffect(){
	
	for (int row=0; row<HEIGHT; row++) {
      for (int col=0; col<WIDTH; col++) {
		  output_buffer[row][col].r = (output_buffer[row][col].r + output_buffer[row][col].g + output_buffer[row][col].b)/3;
		  output_buffer[row][col].g = (output_buffer[row][col].r + output_buffer[row][col].g + output_buffer[row][col].b)/3;
		  output_buffer[row][col].b = (output_buffer[row][col].r + output_buffer[row][col].g + output_buffer[row][col].b)/3;
		  output_buffer[row][col].a = 1.0;
	  }
  }
}

//Function to get user input for selecting the effect,handle the effects and call them 
void Effectsfunc(){
	
	cout<<"You have entered the effects mode, choose any one of the effect to be applied"<<endl;
	cout<<"Now select any of these 6 effects: 1.Overlay effect 2.Warm & cool effect 3. Black n White 4.Deserted effect 5.mirror effect 6.Pray for paris "<<endl;
	int effectchoice;
	scanf("%d",&effectchoice);
	
	if(effectchoice == 1){
		overlay();
	}
	else if(effectchoice == 2){
		wceffects();
	}
	else if(effectchoice == 3){
		greyscaleeffect();
	}
	else if(effectchoice == 4){
		sunshine();
	}
	else if(effectchoice == 5){
		mirroreffect();
	}
	else if(effectchoice == 6){
		paris();
	}
	else {
	cout<<"Choose anyone from 1-6 "<<endl;
	}
	drawImage();
}

//Function to handle perform RGB channel mixing ,based on user input and UP/DOWN keys
void channelmix(int chn,int direction){	
	
	if(direction==1){
		r=1;
		g=1;
		b=1;
		if(chn==1){
			r=r+0.05;
			redflag=1;
			greenflag=0;
			blueflag=0;
		}
		else if(chn==2){
			g=g+0.05;
			redflag=0;
			greenflag=1;
			blueflag=0;
		}
		else if(chn==3){
			b=b+0.05;
			redflag=0;
			greenflag=0;
			blueflag=1;
		}
	}
	else if(direction==0){
		r=1;
		g=1;
		b=1;
		if(chn==1){
			r=r-0.05;
			redflag=1;
			greenflag=0;
			blueflag=0;
		}
		else if(chn==2){
			g=g-0.05;
			redflag=0;
			greenflag=1;
			blueflag=0;
		}
		else if(chn==3){
			b=b-0.05;
			redflag=0;
			greenflag=0;
			blueflag=1;
		}
	}
	else{
		cout<<"wrong direction"<<endl;
	}
	cout<<r<<endl;
	 for(int row = 0; row < HEIGHT; row++){
		for(int col = 0; col < WIDTH; col++) {
			output_buffer[row][col].r = r*(output_buffer[row][col].r);
			output_buffer[row][col].g = g*(output_buffer[row][col].g);
			output_buffer[row][col].b = b*(output_buffer[row][col].b);
			output_buffer[row][col].a = 1.0;
			}
		}  	                                      
	drawImage();
}

//Function to get the user input for selecting the RGB channel to adjust
void channels(){
	
	int channelchoice;
	cout<<"Select the channel you want to increase or decrease 1.RED 2.GREEN 3.BLUE"<<endl;
	scanf("%d",&channelchoice);
	
	if(channelchoice==1)
	{
		redflag=1;
		greenflag=0;
		blueflag=0;
		cout<<" red channel can be altered"<<endl;
	}
	else if(channelchoice==2)
	{
		redflag=0;
		greenflag=1;
		blueflag=0;
		cout<<" green channel can be altered"<<endl;
	}
	else if(channelchoice==3)
	{
		redflag=0;
		greenflag=0;
		blueflag=1;
			cout<<" blue channel can be altered"<<endl;

	}		
	cout<<"Now use up and down arrow keys to increase or decrease"<<endl;
}

//Function to clear all flags before calling a new operation/function
void clearflags(){
	
	gammaflag = 0;
    redflag = 0;
    greenflag = 0;
    blueflag = 0;
    hueflag = 0;
    saturationflag = 0;
    valueflag = 0;
}

//Function to handle UP & DOWN arrow key press
void keycallfunc(int key,int x, int y){
	
	switch (key) {
		case GLUT_KEY_UP :
//Increases the gamma value on up key press							
							if ( gammaflag == 1){
								gammacorrection(1);
								drawImage();
						    }
//Calls the Red,Green or Blue increase based on the choice selected 						    
						    else if (redflag == 1){
								channelmix(1,1);
								cout<<"red call"<<endl;
								drawImage();
						    }
							else if (greenflag == 1){
								channelmix(2,1);
								cout<<"green call"<<endl;
								drawImage();
						    }
						    else if (blueflag == 1){
								channelmix(3,1);
								cout<<"blue call"<<endl;
								drawImage();
						    }
//Calls the Hue,Saturation and Value increase based on the choice selected 					    
						    else if (hueflag == 1){
								HSVcompute(1,1);
								cout<<"hue call"<<endl;
								drawImage();
						    }
							else if (saturationflag == 1){
								HSVcompute(2,1);
								cout<<"saturation call"<<endl;
								drawImage();
						    }
						    else if (valueflag == 1){
								HSVcompute(3,1);
								cout<<"value call"<<endl;
								drawImage();
						    }
						    else{
								cout<<"invalid call"<<endl;
							}						
							break;
		
		case GLUT_KEY_DOWN :							 
//Decreases the gamma value on down key press				
							if ( gammaflag == 1){
								gammacorrection(0);
								drawImage();
						    }
//Calls the Red,Green or Blue increase based on the choice selected 						    
						    else if (redflag == 1){
								channelmix(1,0);
								cout<<"red call"<<endl;
								drawImage();
						    }
							else if (greenflag == 1){
								channelmix(2,0);
								cout<<"green call"<<endl;
								drawImage();
						    }
						    else if (blueflag == 1){
								channelmix(3,0);
								cout<<"blue call"<<endl;
								drawImage();
						    }
//Calls the Hue,Saturation and Value increase based on the choice selected 						    
						    else if (hueflag == 1){
								HSVcompute(1,0);
								cout<<"hue call"<<endl;
								drawImage();
						    }
							else if (saturationflag == 1){
								HSVcompute(2,0);
								cout<<"saturation call"<<endl;
								drawImage();
						    }
						    else if (valueflag == 1){
								HSVcompute(3,0);
								cout<<"value call"<<endl;
								drawImage();
						    }
							else
							{
								cout<<"invalid call"<<endl;
							}
							break;
					}
}	


//Function to handle keypress functions 
void handleKey(unsigned char key, int x, int y){
	
M.identity();
	switch (key)
	{				
//To write the displayed image into the output file    										
	case 'w':	writeimg(Write_file); 			
				break;  
	case 'p':   cout<<"PRESET MODE : Select-"<<endl; 
				cout<<" 'e' for effects"<<endl;
				cout<<" 'f' for setting up frames"<<endl;	 
				break;
	case 'a':   cout<<"ADJUSTMENT MODE : Select-"<<endl;	
				cout<<" 'o' for changing Orientation"<<endl;	
				cout<<" 'h' for adjusting Hue,Saturation and value of the image"<<endl;
				cout<<" 'g' for adjusting the gamma value of the image"<<endl;		
				cout<<" 'c' for adjusting the RGB channels of the image"<<endl;	
				break;
	case 'e': 	backup();			  
				Effectsfunc();
				cout<<"do you want to add another effect? if yes click the image and press e"<<endl;
				break;	
	case 'f': 	cout<<" Add Frame selected "<<endl;
				backup();
				setframes();
				drawImage();
				break;
	case 'o': 	cout<<"Orientation adjust selected"<<endl;
				backup();
				cout<<" Select the type of orientation 1.Horizontal flip 2.Vertical flip 3.Rotate "<<endl;
				scanf("%d", &orientationchoice);
				Orientation(orientationchoice);
				drawImage();
				break;
	case 'h': 	cout<<"Hue,Sat,Value Adjust mode"<<endl;
				backup();
				clearflags();
				HSV();
				break;
	case 'g': 	cout<<"Gamma correction selected"<<endl;
				backup();
				clearflags();
				gammacorrection(3);
				break;
	case 'c': 	cout<<"RGB channels mixer selected"<<endl;
				backup();
				clearflags();
				channels();
				break;
	 case 'u': 	cout<<"Undo option selected"<<endl;
				restorebackup();
				drawImage();
				break;
	case 'x' : 	cout<<"Comparing image with original"<<endl;
				if(compareflag==1){
					drawOriginal();
					compareflag=0;
					}
				else if(compareflag==0){
					drawImage();
					compareflag=1;
					}
		 		break;	            
	case 'q':
	case 'Q':					            
	case 27:  
			   exit(0);
	}
cout<<" Select 'p' for PRESET MODE"<<endl;
cout<<" Select 'a' for ADJUSTMENT MODE"<<endl;
}

// Function to initialize raw data, initiate the load, and call the 
void init_data(){              

//Function call to load image
	loadImage(INPUT_FILE); 
	    
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	frame_buffer = new rgba_pixel*[HEIGHT];
	frame_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      frame_buffer[i] = frame_buffer[i-1] + WIDTH;
	}
	undo_buffer = new rgba_pixel*[HEIGHT];
	undo_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      undo_buffer[i] = undo_buffer[i-1] + WIDTH;
	}  
	flare_buffer = new rgba_pixel*[HEIGHT];
	flare_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	for (int i=1; i<HEIGHT; i++) {
      flare_buffer[i] = flare_buffer[i-1] + WIDTH;
	}  
	output_buffer = new rgba_pixel*[HEIGHT];
    output_buffer[0] = new rgba_pixel[WIDTH*HEIGHT];
	   for (int i=1; i<HEIGHT; i++) {
		  output_buffer[i] = output_buffer[i-1] + WIDTH;
	   }
//function call to populate the output buffer
	populate_output_buffer();           
 
}

//	T.H.E      M.A.I.N        F.U.N.C.T.I.O.N
int main(int argc, char** argv)
{                     
   if (argc < 2) {
      cerr << "usage: " << argv[1] << " [input image]" << endl;
      return -1;
   }
//Copying the write file if it is present    
   strcpy(INPUT_FILE,argv[1]);
   if(argc == 3)
	{
		strcpy(Write_file,argv[2]);     
    }
               
	cout<<"					***** PHOTO EDITOR *****"<<endl;
	cout<<"------------------------------------------------------------------------------"<<endl;
	cout<<"PRESET MODE 'p' : Select-"<<endl; 
	cout<<"			'e' for effects"<<endl;
	cout<<" 			'f' for setting up frames"<<endl;	 
	
	cout<<"ADJUSTMENT MODE 'a': Select-"<<endl;	
	cout<<" 			'o' for changing Orientation"<<endl;	
	cout<<" 			'h' for adjusting Hue,Saturation and value of the image"<<endl;
	cout<<" 			'g' for adjusting the gamma value of the image"<<endl;		
	cout<<" 			'c' for adjusting the RGB channels of the image"<<endl;	
	cout<<" -----------------------------------------------------------------------------"<<endl;	
	cout<<" 'u' for undo the previous action on the image"<<endl;	
	cout<<" 'x' for comparing the image with the original one"<<endl;
	cout<<""<<endl;
	cout<<"***Click on the image and press the keys specified***"<<endl;
	cout<<""<<endl;	

//start up the glut utilities  
   glutInit(&argc, argv);
//Function call to initialize the data
	init_data();  
	drawImage();
//create the graphics window, giving width, height, and title text                  							 
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);  		
	glutInitWindowSize(WIDTH, HEIGHT);	
	glutCreateWindow("PHOTO EDITOR");
//display callback
	glutDisplayFunc(drawImage);	  						
//keyboard callback
	glutKeyboardFunc(handleKey);
	glutSpecialFunc(keycallfunc);	  						
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, WIDTH, 0, HEIGHT);
	glClearColor(1,1,1,1);          
	glutMainLoop();
	return 0;
}


