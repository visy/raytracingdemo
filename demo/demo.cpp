#include <WinSock2.h>
#include <Ws2tcpip.h>
#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "aarect.h"
#include "bvh.h"
#include "pdf.h"
#include "my_stb_obj_loader.h"
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
int max_depth = 5;
float aperture = 0.0;
float dist_to_focus = 100.0;
int samples_per_pixel = 4;
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

// 	spheres.add(make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_center));

	float ss = 1;
	float xxo = -4;
	float yyo = -4;
	float zzo = 8;
    for (int zz = 15; zz > 3; zz--) {
	    int zo = 9-zz;
	    if (zz > 9) {
		 zo = zz-9; 
	    }
	    for (int yy = zo; yy < 9-zo; yy++) {
		    for (int xx = zo; xx < 9-zo; xx++) {
				 float re = cos((zz+yy)*0.4)*0.1;
				 float ge = sin((xx+zz)*0.2)*0.1;
				 float be = cos((xx+yy)*0.4)*0.1;
				 auto fuzz = random_double(0, 0.0);
			auto material_right  = make_shared<metal>(color(1.0-re,1.0-ge,1.0-be), fuzz);
			spheres.add(make_shared<sphere>(point3( xxo+xx*(ss),    zzo-zz*ss, yyo+yy*(ss)),   0.5, material_right));
		    }
	    }
    }
      
    hittable_list ff;
    ff.add(make_shared<bvh_node>(spheres, 0, 1));

    return ff;
}

hittable_list final_scene_lights(){
    // an empty hittable list gives importance sampling with cosine distribution over the normal-hemisphere
    hittable_list lights;
    return lights;
}

int bg_prev = 0;

#define SERVER "valot.instanssi"
#define BUFLEN 512  // max length of answer
#define PORT 9909  // the port on which to listen for incoming data

int vr[22] = {0};
int vg[22] = {0};
int vb[22] = {0};

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

	stream = BASS_StreamCreateFile(FALSE, "music.ogg", 0, 0, BASS_STREAM_PRESCAN);
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
	const sync_track *bg = sync_get_track(rocket, "bg");

	Pixie::Window window;

	int image_height = static_cast<int>(image_width / aspect_ratio);

	w = image_width;
	h = image_height;

	prev = (color*)malloc(w*h*sizeof(color));

	raytraced = (uint32_t*)malloc(w*h*sizeof(uint32_t));
	memset(raytraced,0,w*h);

	// world
	//
	world=final_scene();
	lights = make_shared<hittable_list>(final_scene_lights());

	auto background_sky = make_shared<image_texture>("sky.hdr");
	auto background_moon = make_shared<image_texture>("moon.hdr");
	auto background_street = make_shared<image_texture>("street.hdr");
	background = background_sky;
	background_pdf = make_shared<image_pdf>(background_sky);
	

	// cam
	

	if (!window.Open("Demo or die", xres, yres, DEMO_FULLSCREEN)) 
	{
		return 0;
	}

	WSADATA wsaData;
	int res = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (res != NO_ERROR) {
	    return 0;
	}

	sockaddr_in server;
	int client_socket;
	if ((client_socket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) == SOCKET_ERROR) // <<< UDP socket
	{
		return 2;
	}

	// setup address structure
	memset((char*)&server, 0, sizeof(server));
	server.sin_family = AF_INET;
	server.sin_port = htons(PORT);
	server.sin_addr.S_un.S_addr = inet_addr(SERVER);


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
		}
		if (fe == 1) 
		{
			xo = 1;
		}

		fe++;

		fe = fe % 2;

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
		
		int bg_value =(int(sync_get_val(bg, row)));
		
		lookfrom.x(f_x);
		lookfrom.y(f_y);
		lookfrom.z(f_z);
		
		lookat.x(a_x);
		lookat.y(a_y);
		lookat.z(a_z);

		cam.set_transform(lookfrom,lookat,float(sync_get_val(cam_fov, row)));

		if (bg_value == 0 && bg_prev != bg_value) {
			background = background_sky;
			background_pdf = make_shared<image_pdf>(background_sky);
			bg_prev = bg_value;
		}
		if (bg_value == 1 && bg_prev != bg_value) {
			background = background_moon;
			background_pdf = make_shared<image_pdf>(background_moon);
			bg_prev = bg_value;
		}
		if (bg_value == 2 && bg_prev != bg_value) {
			background = background_street;
			background_pdf = make_shared<image_pdf>(background_street);
			bg_prev = bg_value;
		}
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

		#pragma omp parallel for schedule(dynamic)
		for(int y=0;y<202;y+=1) 
		{
			for(int x=0;x<360;x+=1)
			{
				int i = y*360+x;

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


		int scs[44] = {
			140,20,
			100,20,
			20,20,
			0,0,
			0,40,
			0,80,
			0,140,
			0,180,
			20,180,
			40,180,
			60,180,
			100,180,
			140,180,
			220,180,
			280,180,
			280,140,
			280,80,
			280,40,
			280,20,
			220,20,
			140,20,
			100,20
		};

		for(int j=0;j<22;j++) {
			int rinc = 0;
			int ginc = 0;
			int binc = 0;
			int zx = scs[(j*2)];
			int zy = scs[(j*2)+1];
			for(int yy=0;yy<8;yy++) {
			for(int xx=0;xx<8;xx++) {
			int i = (zy+yy)*360+(zx+xx);
			color pixel = prev[i];
			rinc += clamp(pow(pixel.x()*scale, inv_gamma),0.0,0.999)*256; 
			ginc += clamp(pow(pixel.y()*scale, inv_gamma),0.0,0.999)*256; 
			binc += clamp(pow(pixel.z()*scale, inv_gamma),0.0,0.999)*256; 
			}
			}
			vr[j] = rinc/64;
			vg[j] = ginc/64;
			vb[j] = binc/64;
		}
		bilinear(pixels,raytraced,w,h,xres,yres);

		// update lights
		unsigned char message[157] = { 
		1, 0, 113, 117, 97, 100, 0, 
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0,
		0,0,0,0,0,0};

		for(int i = 0;i<22;i++) 
		{
			message[7+(i*6)+0] = 1; // valo
			message[7+(i*6)+1] = i; // index
			message[7+(i*6)+2] = 0; // extension
			message[7+(i*6)+3] = vr[i]; // r
			message[7+(i*6)+4] = vg[i]; // g
			message[7+(i*6)+5] = vb[i]; // b

		}
		
		sendto(client_socket, message, 157, 0, (sockaddr*)&server, sizeof(sockaddr_in));
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
	closesocket(client_socket);
	WSACleanup();
	window.Close();
}
