/*
CSCI 420
Assignment 3 Raytracer

Name: Rajasimha Reddy Jerry Sivaram Reddy
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pic.h"

// For Linux
// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glut.h>

/* For Mac */
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 200

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 90.0

//defining constants
#define PI 3.14

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//defining a struct point for the vectors
struct point {
   double x;
   double y;
   double z;
};

//vector operations
point normalize(point p){
	double length = sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
	point r;
	r.x = p.x / length;
	r.y = p.y / length;
	r.z = p.z / length;
	return r;
}
point crossProduct(point a, point b){
	point r;
	r.x = a.y * b.z - a.z * b.y;
	r.y = a.z * b.x - a.x * b.z;
	r.z = a.x * b.y - a.y * b.x;
	return r;
}

point addVectors(point a, point b){
	point r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	r.z = a.z + b.z;
	return r;
}

point scalarMultiply(double s, point a){
	point r;
	r.x = a.x * s;
	r.y = a.y * s;
	r.z = a.z * s;
	return r;
}

double dotProduct(point a,point b) {
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

point vectorSubstract(point a, point b) {
  point r;
  r.x = a.x - b.x;
	r.y = a.y - b.y;
	r.z = a.z - b.z;
	return r;
}
//vector operation ends

//A ray structure
struct Ray
{
  point origin;
  point direction;

  bool sphereIntersection(Sphere &sphere, point &intersect);
  bool traingleIntersection(Triangle &triangle , point &intersect);
} ;

//Function to check sphere intersection
bool Ray::sphereIntersection(Sphere &sphere, point &intersect){

  point distance = {origin.x-sphere.position[0],origin.y-sphere.position[1],origin.z-sphere.position[2]};

  double t0,t1;

  double a = dotProduct(direction,direction);
  double b = 2*dotProduct(direction,distance);
  double c = dotProduct(distance,distance) - pow(sphere.radius,2);
  double disc = pow(b,2) - (4*a*c);
  if(disc<0) return false;

  else if (abs(disc) < 1e-8) {
    double q = -0.5f * b/a;
    t0 = q;
    t1 = q;
  }

  else {
    double q = (b>0) ? -0.5f * (b+sqrt(disc)) : -0.5*(b-sqrt(disc));
    t0 = q/a;
    t1 = c/q;
  }

  if(t0<0 && t1<0) return false;
  if(t0>t1 && t1>0) t0=t1;

  intersect = addVectors(origin,scalarMultiply((float)t0,direction));
  
  return true;
}

//Function to check traingle intersection
bool Ray::traingleIntersection(Triangle &triangle , point &intersect) {
  point a = {triangle.v[0].position[0],triangle.v[0].position[1],triangle.v[0].position[2]};
  point b = {triangle.v[1].position[0],triangle.v[1].position[1],triangle.v[1].position[2]};
  point c = {triangle.v[2].position[0],triangle.v[2].position[1],triangle.v[2].position[2]};

  point normal = normalize(crossProduct(vectorSubstract(b,a),vectorSubstract(c,a)));
  double angle = dotProduct(normal,direction);

  if(abs(angle) < 1e-8) return false;

  point distance = vectorSubstract(a,origin);
  double dist = dotProduct(distance,normal);

  double t = dist/angle;
  if(t<0) return false;

  intersect = addVectors(origin,scalarMultiply((float)t,direction));

  if(dotProduct(normal,crossProduct(vectorSubstract(b,a),vectorSubstract(intersect,a))) < 0 ||
  dotProduct(normal,crossProduct(vectorSubstract(c,b),vectorSubstract(intersect,b))) < 0  ||
  dotProduct(normal,crossProduct(vectorSubstract(a,c),vectorSubstract(intersect,c))) < 0) return false;

  return true;
}

void clamp(double &num,float low, float high) {
  if (num>high) num=high;
  else if (num<low) num=low;
}

//Returns the phong shading for the sphere
point getLightingSphere(Sphere &sphere, Light &light, point &intersect){
  point direc1 = {intersect.x-sphere.position[0],intersect.y-sphere.position[1],intersect.z-sphere.position[2]};
  point normal  = normalize(direc1);

  point direc2 = {light.position[0] - intersect.x,light.position[1] - intersect.y,light.position[2] - intersect.z};
  point light_direction = normalize(direc2);

  double magnitude = dotProduct(light_direction,normal);
  clamp(magnitude,0,1);
  point normal_intersect = normalize(scalarMultiply(-1.0f,intersect));
  point reflect = {2*magnitude*normal.x-light_direction.x,2*magnitude*normal.y-light_direction.y,2*magnitude*normal.z-light_direction.z};

  point ref_normal = normalize(reflect);
  double ref_magnitude = dotProduct(ref_normal,normal_intersect);
  clamp(ref_magnitude,0.0f,1.0f);

  point diffuse = {sphere.color_diffuse[0],sphere.color_diffuse[1],sphere.color_diffuse[2]};
  point spec = {sphere.color_specular[0],sphere.color_specular[1],sphere.color_specular[2]};

  point phong_shading = {light.color[0] * (diffuse.x * magnitude + (spec.x * pow(ref_magnitude,sphere.shininess))),
  light.color[1] * (diffuse.y * magnitude + (spec.y * pow(ref_magnitude,sphere.shininess))),
  light.color[2] * (diffuse.z * magnitude + (spec.z * pow(ref_magnitude,sphere.shininess)))};

  return phong_shading;
}

//Returns the phong shading of the triangle
point getLightingTriangle(Triangle &triangle , Light &light, point &intersect) {
  point direc1 = {light.position[0]-intersect.x,light.position[1]-intersect.y,light.position[2]-intersect.z};
  point light_direction = normalize(direc1);

  point a = {triangle.v[0].position[0],triangle.v[0].position[1],triangle.v[0].position[2]};
  point b = {triangle.v[1].position[0],triangle.v[1].position[1],triangle.v[1].position[2]};
  point c = {triangle.v[2].position[0],triangle.v[2].position[1],triangle.v[2].position[2]};

  point ab = vectorSubstract(b,a);
  point ac = vectorSubstract(c,a);
  point cb = vectorSubstract(c,b);
  point ca = vectorSubstract(a,c);

  point intersect_dist_b = vectorSubstract(intersect,b);
  point intersect_dist_c = vectorSubstract(intersect,c);

  point cross_cb = crossProduct(cb,intersect_dist_b);
  point cross_ac = crossProduct(ca,intersect_dist_c);

  point plane = crossProduct(ab,ac);
  float denominator = dotProduct(plane,plane);

  double u = dotProduct(plane,cross_cb)/denominator;
  double v = dotProduct(plane,cross_ac)/denominator;
  double w = 1.0-u-v;

  point normal = {
    u*triangle.v[0].normal[0] + v*triangle.v[1].normal[0] + w*triangle.v[2].normal[0],
    u*triangle.v[0].normal[1] + v*triangle.v[1].normal[1] + w*triangle.v[2].normal[1],
    u*triangle.v[0].normal[2] + v*triangle.v[1].normal[2] + w*triangle.v[2].normal[2]
  };

  point norm_normal = normalize(normal);
  double magnitude = dotProduct(light_direction,norm_normal);
  clamp(magnitude,0,1);
  point normal_intersect = normalize(scalarMultiply(-1.0f,intersect));
  point reflect = {2*magnitude*normal.x-light_direction.x,2*magnitude*normal.y-light_direction.y,2*magnitude*normal.z-light_direction.z};

  point ref_normal = normalize(reflect);
  double ref_magnitude = dotProduct(ref_normal,normal_intersect);
  clamp(ref_magnitude,0.0f,1.0f);

  point diffuse = {
    u*triangle.v[0].color_diffuse[0] + v*triangle.v[1].color_diffuse[0] + w*triangle.v[2].color_diffuse[0],
    u*triangle.v[0].color_diffuse[1] + v*triangle.v[1].color_diffuse[1] + w*triangle.v[2].color_diffuse[1],
    u*triangle.v[0].color_diffuse[2] + v*triangle.v[1].color_diffuse[2] + w*triangle.v[2].color_diffuse[2]
  };
  point spec = {
    u*triangle.v[0].color_specular[0] + v*triangle.v[1].color_specular[0] + w*triangle.v[2].color_specular[0],
    u*triangle.v[0].color_specular[1] + v*triangle.v[1].color_specular[1] + w*triangle.v[2].color_specular[1],
    u*triangle.v[0].color_specular[2] + v*triangle.v[1].color_specular[2] + w*triangle.v[2].color_specular[2]
  };
  double shininess = u*triangle.v[0].shininess + v*triangle.v[1].shininess + w*triangle.v[2].shininess;

  point phong_shading = {light.color[0] * (diffuse.x * magnitude + (spec.x * pow(ref_magnitude,shininess))),
  light.color[1] * (diffuse.y * magnitude + (spec.y * pow(ref_magnitude,shininess))),
  light.color[2] * (diffuse.z * magnitude + (spec.z * pow(ref_magnitude,shininess)))};

  return phong_shading;
}

//Function to check collision
point getCollision(Ray &ray, point &color, double &near) {
  point finColor = color;

  for(int i=0;i<num_spheres;i++) {
    point intersect={0,0,-1e8};
    bool isIntersecting=ray.sphereIntersection(spheres[i],intersect);
    if(isIntersecting&&intersect.z>near) {
      finColor.x = finColor.y = finColor.z =  0.0f;
      for(int j=0;j<num_lights;j++) {
        bool isLit = true;
        point lightPosition = {lights[j].position[0],lights[j].position[1],lights[j].position[2]};
        Ray shadow;
        shadow.origin = intersect;
        shadow.direction = normalize(vectorSubstract(lightPosition,intersect));
        for(int k=0;k<num_spheres;k++){
          point shadowIntersect;
          if(shadow.sphereIntersection(spheres[k],shadowIntersect) && k!=i) {
            point a = vectorSubstract(shadowIntersect,intersect);
            point b = vectorSubstract(lightPosition,intersect);
            if(sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) < sqrt(pow(b.x, 2) + pow(b.y, 2) + pow(b.z, 2))) {
              isLit = false;
              break;
            }
          }
        }

        for(int k=0;k<num_triangles;k++){
          point shadowIntersect;
          if(shadow.traingleIntersection(triangles[k],shadowIntersect)) {
            point a = vectorSubstract(shadowIntersect,intersect);
            point b = vectorSubstract(lightPosition,intersect);
            if(sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) < sqrt(pow(b.x, 2) + pow(b.y, 2) + pow(b.z, 2))) {
              isLit = false;
              break;
            }
          }
        }

      if(isLit) finColor = addVectors(finColor,getLightingSphere(spheres[i],lights[j],intersect));

      }
      near = intersect.z;
    }
  }

  for(int i=0;i<num_triangles;i++) {
    point intersect={0,0,-1e8};
    bool isIntersecting=ray.traingleIntersection(triangles[i],intersect);
    if(isIntersecting&&intersect.z>near) {
      finColor.x = finColor.y = finColor.z =  0.0f;
      for(int j=0;j<num_lights;j++) {
        bool isLit = true;
        point lightPosition = {lights[j].position[0],lights[j].position[1],lights[j].position[2]};
        Ray shadow;
        shadow.origin = intersect;
        shadow.direction = normalize(vectorSubstract(lightPosition,intersect));
        for(int k=0;k<num_spheres;k++){
          point shadowIntersect;
          if(shadow.sphereIntersection(spheres[k],shadowIntersect)) {
            point a = vectorSubstract(shadowIntersect,intersect);
            point b = vectorSubstract(lightPosition,intersect);
            if(sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) < sqrt(pow(b.x, 2) + pow(b.y, 2) + pow(b.z, 2))) {
              isLit = false;
              break;
            }
          }
        }

        for(int k=0;k<num_triangles;k++){
          point shadowIntersect;
          if(shadow.traingleIntersection(triangles[k],shadowIntersect)  && k!=i) {
            point a = vectorSubstract(shadowIntersect,intersect);
            point b = vectorSubstract(lightPosition,intersect);
            if(sqrt(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2)) < sqrt(pow(b.x, 2) + pow(b.y, 2) + pow(b.z, 2))) {
              isLit = false;
              break;
            }
          }
        }

        if(isLit) finColor = addVectors(finColor,getLightingTriangle(triangles[i],lights[j],intersect));
      }
      near = intersect.z;
    }
  }

  return finColor;
}

point reflect(point &I,point &N) {
  return(vectorSubstract(I,scalarMultiply(2,scalarMultiply(dotProduct(I,N),N))));
}

//Returns the color of each ray
point getRay(Ray &ray){
  point color={1.0f,1.0f,1.0f};
  int i=4;
  double intersect = -1e8;
  color = scalarMultiply(0.8f,getCollision(ray,color,intersect));
  //if(color.x == 1&&color.y==1&&color.z==1) return color;
  // while(i>=0) {
  //   point intersectPoint = {0,0,intersect};
  //   point intersectNormal = normalize(intersectPoint);
  //   point R = reflect(ray.direction,intersectNormal);
  //   ray.origin = addVectors(intersectPoint,intersectNormal);
  //   ray.direction = R;
  //   color = addVectors(color,scalarMultiply(0.8f,getCollision(ray,color,intersect)));
  //   i--;
  // }
  //color = getCollision(ray,color,intersect);
  //printf("intersect=%f\n",intersect);
  point ambient = {ambient_light[0],ambient_light[1],ambient_light[2]};
  color = addVectors(color,ambient);
  return color;
}

//Get the camera ray of each point
Ray getCameraRay(double x , double y){
  double aspectRatio = (double)WIDTH/(double)HEIGHT;
  double angle = tan((fov/2)*PI/180);
  double xCamera = ((2*(x+0.5)/(double)WIDTH)-1)*aspectRatio*angle;
  double yCamera = ((2*(y+0.5)/(double)HEIGHT)-1)*aspectRatio;
  point origin = {0.0f,0.0f,0.0f};
  point direction = {xCamera,yCamera,-1};
  Ray ray;
  ray.origin = origin;
  ray.direction = normalize(direction);
  return ray;
}

//get random light rays for the soft shadow
void randomLights(int  m, int finalLight, double num_lights)
{
	double delta = 0.05;

	static int randomNumbers[] = { -4, -3, -2, -1, 1, 2, 3, 4} ;
	int index = rand() % (sizeof randomNumbers / sizeof *randomNumbers);
	int random = randomNumbers[index];

	lights[finalLight].position[0] = lights[m].position[0] + random * delta;
	lights[finalLight].position[1] = lights[m].position[1] + random * delta;
	lights[finalLight].position[2] = lights[m].position[2] + random * delta;

	lights[finalLight].color[0] = lights[m].color[0] / num_lights;
	lights[finalLight].color[1] = lights[m].color[1] / num_lights;
	lights[finalLight].color[2] = lights[m].color[2] / num_lights;
}

//add those lights to the light array
void addLights()
{
	int finalLight = num_lights;

	for (int m = 0; m < num_lights; m++)
	{
		for (int i = 0; i < 32; i++)
		{
			randomLights(m, finalLight++, 32);
		}
		lights[m].color[0] = lights[m].color[0] / 32;
		lights[m].color[1] = lights[m].color[1] / 32;
		lights[m].color[2] = lights[m].color[2] / 32;
	}

	num_lights = finalLight;
}

//MODIFY THIS FUNCTION

//Driver code
void draw_scene()
{
  bool antialiasing = false;
  bool softShadows = true;
  unsigned int x,y;
  //simple output

  //for soft shadows
  if(softShadows) addLights();

  for(x=0; x<WIDTH; x++)
  {
    glPointSize(1.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      point color;
      color.x = color.y = color.z = 0.0f;

      //Antialiasing
      if(antialiasing) {
        for (int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
          Ray ray = getCameraRay(x+i,y+j);
          point rayColor = getRay(ray);
          color = addVectors(color,rayColor);
        }
      }
      color = scalarMultiply(0.1111f,color);
      }
      //No Antialiasing
      else {
      Ray ray = getCameraRay(x,y);
      color = getRay(ray);
      }
      plot_pixel(x,y,color.x*255,color.y*255,color.z*255);

    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
