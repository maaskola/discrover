#include <iostream>
#include <cairo.h>
#include "logo.hpp"

using namespace std;

bool draw_logo(const string &path) {
  string out_path = "hello.png";
  cout << "Drawing the logo for " << path << " to " << out_path << "." << endl;

  cairo_surface_t *surface =
    cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 240, 80);
  cairo_t *cr =
    cairo_create (surface);

  cairo_select_font_face (cr, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size (cr, 32.0);
  cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);
  cairo_move_to (cr, 10.0, 50.0);
  cairo_show_text (cr, "Hello, world");

  cairo_destroy (cr);
  cairo_surface_write_to_png (surface, out_path.c_str());
  cairo_surface_destroy (surface);
  return EXIT_SUCCESS;
}
