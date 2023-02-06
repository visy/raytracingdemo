#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "bvh.h"
#include "pdf.h"
#include "pixie.h"
#include "bass.h"
#include "sync.h"

#include <omp.h>
#include <time.h>

static const float bpm = 150.0f; /* beats per minute */
static const int rpb = 8; /* rows per beat */
static const double row_rate = (double(bpm) / 60) * rpb;

static double bass_get_row(HSTREAM h)
{
	QWORD pos = BASS_ChannelGetPosition(h, BASS_POS_BYTE);
	double time = BASS_ChannelBytes2Seconds(h, pos);
	return time * row_rate;
}

#define SYNC_PLAYER

#ifndef SYNC_PLAYER

static void bass_pause(void *d, int flag)
{
	HSTREAM h = *((HSTREAM *)d);
	if (flag)
		BASS_ChannelPause(h);
	else
		BASS_ChannelPlay(h, false);
}

static void bass_set_row(void *d, int row)
{
	HSTREAM h = *((HSTREAM *)d);
	QWORD pos = BASS_ChannelSeconds2Bytes(h, row / row_rate);
	BASS_ChannelSetPosition(h, pos, BASS_POS_BYTE);
}

static int bass_is_playing(void *d)
{
	HSTREAM h = *((HSTREAM *)d);
	return BASS_ChannelIsActive(h) == BASS_ACTIVE_PLAYING;
}

static struct sync_cb bass_cb = {
	bass_pause,
	bass_set_row,
	bass_is_playing
};

#endif



int w = 0;
int h = 0;
float t = 0;
float c_r,c_g,c_b=1;

color ray_color(const ray& r, const shared_ptr<texture>& background, const shared_ptr<pdf>& background_pdf, const hittable& world, const shared_ptr<hittable>& lights, int depth) 
{
	hit_record rec;

	if (depth <= 0) 
	{
		return color(0,0,0);
	}

	if(!world.hit(r, 0.000001, infinity, rec)){
		auto unit_dir = unit_vector(r.direction());
		double u, v; get_spherical_uv(unit_dir, u, v);
		return background->value(u, v, unit_dir);
	}

	scatter_record srec;
	color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);

	if (!rec.mat_ptr->scatter(r, rec, srec)){
		return emitted;
	}

	// no importance sampling
	if (srec.is_specular){
		return srec.attenuation * ray_color(srec.specular_ray, background, background_pdf, world, lights, depth-1);
	}

	auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
	mixture_pdf p_objs(light_ptr, srec.pdf_ptr, 0.5);
	mixture_pdf p(make_shared<mixture_pdf>(p_objs), background_pdf, 0.8);

	ray scattered = ray(rec.p, p.generate(), r.time());
	auto pdf_val = p.value(scattered.direction());

	return emitted + 
	srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
			 * ray_color(scattered, background, background_pdf, world, lights, depth-1) / pdf_val;
}

void bilinear(uint32_t* pd, uint32_t* pixels, int widthSource, int heightSource, int width, int height)
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
                r = (byte)((rf));
                g = (byte)((gf));
                b = (byte)((bf));

                pd[srcIdx++] = MAKE_RGB(r,g,b);
        }
    }
}

color* prev;
uint32_t* raytraced;

int fe,xo,yo = 0;
hittable_list world;
shared_ptr<hittable_list> lights;
point3 lookfrom(-2,2,5);
point3 lookat(0,0,-1);
int xres = 1280;
int yres = 720;
auto aspect_ratio = 16.0 / 9.0;
int image_width = 360;
int max_depth = 6;
float aperture = 0.0;
float dist_to_focus = 10.0;
int samples_per_pixel = 8;

camera cam(lookfrom, lookat, vec3(0,1,0), 75, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

shared_ptr<texture> background;
shared_ptr<pdf> background_pdf;	

color render_samples(int x, int y, int w, int h) {
	color pixel(0,0,0);
	for (int i = 0; i < samples_per_pixel; i++) {
		auto u = (x+random_double()) / (w-1);
		auto v = ((h-y)+random_double()) / (h-1);
		ray r = cam.get_ray(u,v);
		color ray_contrib = ray_color(r, background, background_pdf, world, lights, max_depth);
		zero_nan_vals(ray_contrib);
		pixel += ray_contrib;
	}

	return pixel;
}

hittable_list final_scene() {
	auto material_center = make_shared<lambertian>(color(0.8, 0.3, 0.3));
	auto material_left   = make_shared<dielectric>(1.5);

	auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.2, 0.9));

	hittable_list spheres;

// 	spheres.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, make_shared<lambertian>(checker)));

	float ss = 1;
	float xxo = -4;
	float yyo = -4;
    for (int yy = 0; yy < 8; yy++) {
	    for (int xx = 0; xx < 8; xx++) {
			 float re = cos((xx+yy)*0.5)*0.2;
			 float ge = sin((xx+yy)*0.1)*0.3;
			 float be = cos((xx+yy)*0.2)*0.1;
			 auto fuzz = random_double(0, 0.5);
		auto material_right  = make_shared<metal>(color(0.5-re,0.5-ge,0.5-be), fuzz);
		spheres.add(make_shared<sphere>(point3( xxo+xx*ss,    0.0, yyo+yy*ss),   0.5, material_right));
	    }
    }

        hittable_list objects;
    objects.add(make_shared<bvh_node>(spheres, 0, 1));

    return objects;
}

hittable_list final_scene_lights(){
    // an empty hittable list gives importance sampling with cosine distribution over the normal-hemisphere
    hittable_list lights;
    return lights;
}

int main(int argc, char** argv)
{
	HSTREAM stream;
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

	stream = BASS_StreamCreateFile(FALSE, "music.wav", 0, 0, BASS_STREAM_PRESCAN);
	if (!stream) {
		return 0;
	}

	sync_device *rocket = sync_create_device("sync");

#ifndef SYNC_PLAYER
	if (sync_tcp_connect(rocket, "localhost", SYNC_DEFAULT_PORT)) {
		return 0;
	}
#endif

	const sync_track *clear_r = sync_get_track(rocket, "clear.r");
	const sync_track *clear_g = sync_get_track(rocket, "clear.g");
	const sync_track *clear_b = sync_get_track(rocket, "clear.b");
	const sync_track *cam_from_x = sync_get_track(rocket, "from.x");
	const sync_track *cam_from_y = sync_get_track(rocket, "from.y");
	const sync_track *cam_from_z = sync_get_track(rocket, "from.z");
	const sync_track *cam_at_x = sync_get_track(rocket, "at.x");
	const sync_track *cam_at_y = sync_get_track(rocket, "at.y");
	const sync_track *cam_at_z = sync_get_track(rocket, "at.z");
	const sync_track *cam_fov = sync_get_track(rocket, "fov");
	const sync_track *effu1 = sync_get_track(rocket, "effu1");

	Pixie::Window window;

	int image_height = static_cast<int>(image_width / aspect_ratio);

	w = image_width;
	h = image_height;

	prev = (color*)malloc(w*h*sizeof(color));

	raytraced = (uint32_t*)malloc(w*h*sizeof(uint32_t));

	// world
	//
	world=final_scene();
	lights = make_shared<hittable_list>(final_scene_lights());

	auto background_skybox = make_shared<image_texture>("studio.hdr");
	background = background_skybox;
	background_pdf = make_shared<image_pdf>(background_skybox);
	

	// cam
	

	if (!window.Open("Demo or die", xres, yres, DEMO_FULLSCREEN)) 
	{
		return 0;
	}

	BASS_Start();	
	BASS_ChannelPlay(stream, FALSE);


	while (!window.HasKeyGoneUp(Pixie::Key_Escape))
	{
		double row = bass_get_row(stream);

#ifndef SYNC_PLAYER
		if (sync_update(rocket, (int)floor(row), &bass_cb, (void *)&stream))
			sync_tcp_connect(rocket, "localhost", SYNC_DEFAULT_PORT);
#endif

		uint32_t* pixels = window.GetPixels();
	
		// ..draw pixels!
		float delta = window.GetDelta();
		t=t+delta*1000;
		int tt = (int)t;
				
		float scale = 1.0 / samples_per_pixel;
		auto inv_gamma = 1./2.2f;

		if (fe == 0) 
		{
			xo = 0;
			yo = 0;
		}
		if (fe == 1) 
		{
			xo = 1;
			yo = 1;
		}
		if (fe == 2) 
		{
			xo = 0;
			yo = 1;
		}
		if (fe == 3) 
		{
			xo = 1;
			yo = 0;
		}

		fe++;

		fe = fe % 4;

		// 
		c_r =(float(sync_get_val(clear_r, row)));
		c_g =(float(sync_get_val(clear_g, row)));
		c_b =(float(sync_get_val(clear_b, row)));
		
		float f_x =(float(sync_get_val(cam_from_x, row)));
		float f_y =(float(sync_get_val(cam_from_y, row)));
		float f_z =(float(sync_get_val(cam_from_z, row)));
		
		float a_x =(float(sync_get_val(cam_at_x, row)));
		float a_y =(float(sync_get_val(cam_at_y, row)));
		float a_z =(float(sync_get_val(cam_at_z, row)));
		
		lookfrom.x(f_x);
		lookfrom.y(f_y);
		lookfrom.z(f_z);
		
		lookat.x(a_x);
		lookat.y(a_y);
		lookat.z(a_z);

		cam.set_transform(lookfrom,lookat,float(sync_get_val(cam_fov, row)));
/*		
		float cyc = float(sync_get_val(effu1,row));
		int i = 0;
		 for (auto obu : world.objects) {
			 i++;
			 if (i == 1) continue;
			 float re = cos(i+cyc*0.005)*0.2;
			 float ge = sin(i+cyc*0.001)*0.3;
			 float be = cos(i+cyc*0.002)*0.1;
    			auto material_right  = make_shared<metal>(color(0.7-re, 0.2+ge, 0.8-be));
		 }
*/


		const int N_THREADS = 16;
		const int CHUNKS_PER_THREAD = 3;

		int global_done_scanlines=0;
		#pragma omp parallel num_threads(N_THREADS)
		{
		srand(int(time(NULL))^ omp_get_thread_num());
		#pragma omp for schedule(dynamic, h/(N_THREADS*CHUNKS_PER_THREAD))
		for(int y=yo;y<h;y+=2) 
		{
			
			#pragma omp atomic
			++global_done_scanlines;
			for(int x=xo;x<w;x+=2)
			{
				int i = y*w+x;

				color pixel(0,0,0);

				pixel = render_samples(x,y,w,h);
	
				pixel = (pixel + prev[i])*vec3(0.5,0.5,0.5);
				prev[i] = pixel;
				int r = clamp(pow(pixel.x()*scale, inv_gamma),0.0,0.999)*256; 
				int g = clamp(pow(pixel.y()*scale, inv_gamma),0.0,0.999)*256; 
				int b = clamp(pow(pixel.z()*scale, inv_gamma),0.0,0.999)*256; 
				raytraced[i] = MAKE_RGB(r,g,b);
			}
		}
		}

		bilinear(pixels,raytraced,w,h,xres,yres);

		
		if (!window.Update())
		{
			break;
		}

		BASS_Update(0);
	}
#ifndef SYNC_PLAYER
	sync_save_tracks(rocket);
#endif
	sync_destroy_device(rocket);

	BASS_StreamFree(stream);
	BASS_Free();
	window.Close();
}
