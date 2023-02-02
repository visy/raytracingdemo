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

void bilinear(uint32_t* pd, uint32_t* pixels, int widthSource, int heightSource, int width, int height, byte* baseImage)
{
    float xs = (float)widthSource / width;
    float ys = (float)heightSource / height;

    float fracx, fracy, ifracx, ifracy, sx, sy, l0, l1, rf, gf, bf;
    int c, x0, x1, y0, y1;
    byte c1a, c1r, c1g, c1b, c2a, c2r, c2g, c2b, c3a, c3r, c3g, c3b, c4a, c4r, c4g, c4b;
    byte a, r, g, b;

    // Bilinear
    int srcIdx = 0;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            sx = x * xs;
            sy = y * ys;
            x0 = (int)sx;
            y0 = (int)sy;

            // Calculate coordinates of the 4 interpolation points
            fracx = sx - x0;
            fracy = sy - y0;
            ifracx = 1.0f - fracx;
            ifracy = 1.0f - fracy;
            x1 = x0 + 1;
            if (x1 >= widthSource)
            {
                x1 = x0;
            }
            y1 = y0 + 1;
            if (y1 >= heightSource)
            {
                y1 = y0;
            }

            // Read source color
            c = pixels[y0 * widthSource + x0];
            c1a = (byte)(c >> 24);
            c1r = (byte)(c >> 16);
            c1g = (byte)(c >> 8);
            c1b = (byte)(c);

            c = pixels[y0 * widthSource + x1];
            c2a = (byte)(c >> 24);
            c2r = (byte)(c >> 16);
            c2g = (byte)(c >> 8);
            c2b = (byte)(c);

            c = pixels[y1 * widthSource + x0];
            c3a = (byte)(c >> 24);
            c3r = (byte)(c >> 16);
            c3g = (byte)(c >> 8);
            c3b = (byte)(c);

            c = pixels[y1 * widthSource + x1];
            c4a = (byte)(c >> 24);
            c4r = (byte)(c >> 16);
            c4g = (byte)(c >> 8);
            c4b = (byte)(c);

            // Calculate colors
            // Alpha
            l0 = ifracx * c1a + fracx * c2a;
            l1 = ifracx * c3a + fracx * c4a;
            a = (byte)(ifracy * l0 + fracy * l1);

            // Write destination
            if (a > 0)
            {
                // Red
                l0 = ifracx * c1r + fracx * c2r;
                l1 = ifracx * c3r + fracx * c4r;
                rf = ifracy * l0 + fracy * l1;

                // Green
                l0 = ifracx * c1g + fracx * c2g;
                l1 = ifracx * c3g + fracx * c4g;
                gf = ifracy * l0 + fracy * l1;

                // Blue
                l0 = ifracx * c1b + fracx * c2b;
                l1 = ifracx * c3b + fracx * c4b;
                bf = ifracy * l0 + fracy * l1;

                // Cast to byte
                float alpha = a / 255.0f;
                r = (byte)((rf * alpha) + (baseImage[srcIdx] * (1.0f - alpha)));
                g = (byte)((gf * alpha) + (baseImage[srcIdx] * (1.0f - alpha)));
                b = (byte)((bf * alpha) + (baseImage[srcIdx] * (1.0f - alpha)));

                pd[srcIdx++] = (255 << 24) | (r << 16) | (g << 8) | b;
            }
            else
            {
                // Alpha, Red, Green, Blue                          
                pd[srcIdx++] = (255 << 24) | (baseImage[srcIdx] << 16) | (baseImage[srcIdx] << 8) | baseImage[srcIdx];
            }
        }
    }
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

//	chan = BASS_StreamCreateFile(FALSE, "music.wav", 0, 0, BASS_SAMPLE_FLOAT);
//	if (!chan && BASS_ErrorGetCode() == BASS_ERROR_FILEFORM) {
//		return 0;
//	}


	Pixie::Window window;

	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 512;
	int image_height = static_cast<int>(image_width / aspect_ratio);
	int samples_per_pixel = 16;
	int max_depth = 5;
	int gridsize = 3;

	w = image_width;
	h = image_height;

	prev = (color*)malloc(w*h*sizeof(color));

	// world
	hittable_list world;
	world.add(make_shared<sphere>(point3(0,0,-1), 0.4, color(0.7,0.5,0.0)));
	world.add(make_shared<sphere>(point3(0.75,0,-1), 0.4, color(1.0,0.0,0.0)));
	world.add(make_shared<sphere>(point3(-0.75,0,-1), 0.4, color(0.0,1.0,1.0)));
	world.add(make_shared<sphere>(point3(0,-100.5,-1), 100, color(0.0,0.0,1.0)));

	// cam
	camera cam;
	

	if (!window.Open("Demo or die", w, h, DEMO_FULLSCREEN)) 
	{
		return 0;
	}
	
//	BASS_ChannelPlay(chan, FALSE);


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
//		cam.origin.z(-025+cos(t*0.005)*0.32);
		
		if (!window.Update())
		{
			break;
		}
	}

	BASS_Free();
	window.Close();
}
