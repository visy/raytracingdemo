#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "pixie.h"

#include "bass.h"

int w = 0;
int h = 0;
float t = 0;

color ray_color(const ray& r, const hittable& world, int depth, color c) 
{
	hit_record rec;

	if (depth <= 0) 
	{
		return color(0,0,0);
	}

	if (world.hit(r, 0.001, infinity, rec)) 
	{
		point3 target = rec.p + random_in_hemisphere(rec.normal);
		return 0.5 * ray_color(ray(rec.p, target - rec.p), world, depth-1, rec.c);
	}
	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5*(unit_direction.y() + 1.0);
	return (1.0-t)*color(1.0, 1.0, 1.0) + t*c;
}

color* prev;

int main(int argc, char** argv)
{
	DWORD chan, act, level;
	QWORD pos;
	double secs;
	int a, filep, device = -1;
	BASS_CHANNELINFO info;

	if (HIWORD(BASS_GetVersion()) != BASSVERSION) {
		return 0;
	}

	if (!BASS_Init(device, 48000, 0, 0, NULL)) {
		return 0;
	}

	chan = BASS_StreamCreateFile(FALSE, "music.wav", 0, 0, BASS_SAMPLE_FLOAT);
	if (!chan && BASS_ErrorGetCode() == BASS_ERROR_FILEFORM) {
		return 0;
	}


	Pixie::Window window;

	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 512;
	int image_height = static_cast<int>(image_width / aspect_ratio);
	int samples_per_pixel = 64;
	int max_depth = 5;

	w = image_width;
	h = image_height;

	prev = (color*)malloc(w*h*sizeof(color));

	// world
	hittable_list world;
	world.add(make_shared<sphere>(point3(0,0,-1), 0.3, color(0.7,0.5,0.0)));
	world.add(make_shared<sphere>(point3(0.5,0,-1), 0.3, color(1.0,0.0,0.0)));
	world.add(make_shared<sphere>(point3(-0.5,0,-1), 0.3, color(0.0,1.0,1.0)));
	world.add(make_shared<sphere>(point3(0,-100.5,-1), 100, color(0.0,0.0,1.0)));

	// cam
	camera cam;
	

	if (!window.Open("Demo or die", w, h, DEMO_FULLSCREEN)) 
	{
		return 0;
	}
	
	BASS_ChannelPlay(chan, FALSE);

	int gridsize = 16;

	while (!window.HasKeyGoneUp(Pixie::Key_Escape))
	{
		uint32_t* pixels = window.GetPixels();
		// ..draw pixels!
		float delta = window.GetDelta();
		t=t+delta*1000;
		int tt = (int)t;
				
		float scale = 1.0 / samples_per_pixel;

		#pragma omp for
		for(int y=0;y<h;y+=gridsize) 
		{
			for(int x=0;x<w;x+=gridsize) 
			{
				int i = y*w+x;

				color pixel(0,0,0);

				for (int s = 0;s < samples_per_pixel; s++) 
				{
					auto u = (x+random_double()) / (w-1);
					auto v = ((h-y)+random_double()) / (h-1);
					ray r = cam.get_ray(u,v);
					pixel += ray_color(r, world, max_depth, color(0,0,0));
				}
	
				pixel = (pixel + prev[i])*vec3(0.5,0.5,0.5);
				prev[i] = pixel;
			}
		}

		#pragma omp for
		for(int y=0;y<h-gridsize;y+=gridsize) 
		{
			for(int x=0;x<w-gridsize;x+=gridsize)
			{
				color sample = prev[y*w+x];
				
				color sample2;

				for(int y2=1;y2<gridsize+1;y2+=1)
				{
					for(int x2=1;x2<gridsize+1;x2+=1)
					{
						if (y+gridsize < h) {
							sample2	= prev[(y+y2)*w+(x+x2)];
						}
						int i = ((y+y2)*w)+(x+x2);
						color pixel = lerp(sample,sample2,(x2-1)/((float)gridsize));

						int r = sqrt(pixel.x()*scale)*255; 
						int g = sqrt(pixel.y()*scale)*255; 
						int b = sqrt(pixel.z()*scale)*255; 

						pixels[i] = MAKE_RGB(r,g,b);
					}
				}
			}
		}
		cam.origin.z(-0.25+cos(t*0.005)*0.32);
		
		if (!window.Update())
		{
			break;
		}
	}

	BASS_Free();
	window.Close();
}
