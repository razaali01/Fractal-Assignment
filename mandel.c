/*
	Syed Ali
	1001725463
*/

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>				// for clock_t, clock(), CLOCKS_PER_SEC
#include <unistd.h>	    	 // for sleep()
#include <pthread.h>			// for threads

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void * compute_image(void * arg);


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("-n <threads> Give the amount of threads you would like this run to have");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

// create a structure to hold the values for
// each of the threads
struct parameters
{
	int thread_id;
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	int NumberOfThreads;
};



int main( int argc, char *argv[] )
{
	// variables to hold to start time of the code
	// and the end time of the code
	struct timeval begin_time;
  struct timeval end_time;

	// this is where the code starts so this
	// is where we want to start the clock to
	// keep track of time the code took to
	// execute
	gettimeofday( &begin_time, NULL );


	char c;

	// These are the default configuration values used
	// if no command line arguments are given.
	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int image_width = 500;
	int image_height = 500;
	int max = 1000;
	// Creating a default value of 1 thread if
	// no -n is given
	int NumberOfThreads = 1;

	// For each command line argument given,
	// override the appropriate configuration value.

	// modifacation of command line to include -n
	// while((c = getopt(argc,argv,"x:y:s:W:H:m:o:h"))!=-1)

	while ((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1)
	{
		switch(c)
		{
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
				NumberOfThreads = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Creating the threads and creating the structs
	// based on how many threads are being used
	pthread_t threads[NumberOfThreads];
	struct parameters params[NumberOfThreads];



	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d threads= %d outfile=%s\n",xcenter,ycenter,scale,max,NumberOfThreads,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	// compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max);

	int i;


	// For each thread load the values retireved from the
	// command line into the struct for the thread
	for(i = 0; i < NumberOfThreads; i++)
	{
		params[i].bm = bm;
		params[i].xmin = xcenter-scale;
		params[i].xmax = xcenter+scale;
		params[i].ymin = ycenter-scale;
		params[i].ymax = ycenter+scale;
		params[i].max = max;
		params[i].thread_id = i;
		params[i].NumberOfThreads = NumberOfThreads;

		pthread_create(&threads[i], NULL, compute_image, (void *) &params[i]);
	}

	for(i = 0; i < NumberOfThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}


	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile))
	{
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}


	// this is where the code ends so this
	// is where we want to stop the clock
	gettimeofday( &end_time, NULL );

	// calculate time by difference of end and begin
	// then dividing the diff by CLOCK_PER_SEC to
	// convert to seconds
	long time_to_execute1 = ( end_time.tv_sec * 1000000 + end_time.tv_usec );
	long time_to_execute2 = ( begin_time.tv_sec * 1000000 + begin_time.tv_usec );
	long time_to_execute = time_to_execute1 - time_to_execute2;

	printf("This code took %d microseconds to execute\n", time_to_execute);
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void * compute_image(void * arg)
{
	int i,j;

	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	int NumberOfThreads;


	struct parameters * params = (struct parameters *) arg;
	// get each variables value from the struct
	// params and assign it to variables
	bm = params -> bm;
	xmin = params -> xmin;
	xmax = params -> xmax;
	ymin = params -> ymin;
	ymax = params -> ymax;
	max = params -> max;
	NumberOfThreads = params -> NumberOfThreads;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// calculate which section of the bitmap each of the
	// threads will be working on
	int beginH = params -> thread_id * height / NumberOfThreads;
	int endH = (beginH + height / NumberOfThreads) + 1;

	// For every pixel in the image...

	for(j = beginH;j < endH; j++) {

		for(i = 0; i < width; i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
  return NULL;
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
