#include <iostream>
#include <vector>
#include <cairo.h>
#include <cairo-pdf.h>
#include "logo.hpp"

namespace logo {
using namespace std;

struct coord_t {
  double x;
  double y;
};

struct rgb_t {
  double r;
  double g;
  double b;
};

struct palette_t {
  rgb_t a;
  rgb_t c;
  rgb_t g;
  rgb_t t;
};

palette_t default_colors = {
  {0.0, 0.75, 0.0},
  {0.0, 0.0, 1.0},
  {1, 0.6470588, 0.0},
  {1.0, 0.0, 0.0}
};

using coords_t = vector<coord_t>;

coords_t letter_a = {
  {0, 0},
  {4, 10},
  {6, 10},
  {10, 0},
  {8, 0},
  {6.8, 3},
  {3.2, 3},
  {2, 0},
  {0, 0},
  {3.6, 4},
  {5, 7.5},
  {6.4, 4},
  {3.6, 4}
};

coords_t letter_t = {
  {10, 10},
  {10, 9},
  {6, 9},
  {6, 0},
  {4, 0},
  {4, 9},
  {0, 9},
  {0, 10}
};

void draw_letter(cairo_t *cr, const coords_t &coords, double x, double y, double width, double height, rgb_t color) {
  cairo_move_to(cr, x + coords[0].x / 10 * width, y - coords[0].y / 10 * height);
  for(auto coord: coords)
    cairo_line_to(cr, x + coord.x / 10 * width, y - coord.y / 10 * height);
  cairo_set_line_width (cr, 1.0);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_EVEN_ODD);
  cairo_set_source_rgb (cr, color.r, color.g, color.b);
  cairo_fill(cr);
}

void draw_logo_to_surface(cairo_surface_t *surface) {
  cairo_t *cr =
    cairo_create (surface);

  cairo_select_font_face (cr, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size (cr, 32.0);
  cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);
  cairo_move_to (cr, 10.0, 50.0);
  cairo_show_text (cr, "Hello, world");

  draw_letter(cr, letter_a, 10, 50, 30, 30, default_colors.a);
  draw_letter(cr, letter_t, 40, 50, 30, 30, default_colors.c);
  draw_letter(cr, letter_t, 40, 80, 30, 30, default_colors.g);
  draw_letter(cr, letter_a, 40, 20, 30, 30, default_colors.t);

  cairo_show_page(cr);
  cairo_destroy(cr);
}

string ending(output_t kind) {
  switch(kind) {
    case output_t::PDF:
      return "pdf";
    case output_t::PNG:
      return "png";
    default:
      return "";
  };
}

bool draw_logo(const string &path, output_t kind) {
  string out_path = path + "." + ending(kind);
  cout << "Drawing the logo for " << path << " to " << out_path << "." << endl;

  cairo_surface_t *surface;
  switch(kind) {
    case output_t::PDF:
      surface = cairo_pdf_surface_create (out_path.c_str(), 240, 80);
      break;
    case output_t::PNG:
      surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 240, 80);
      break;
  }

  draw_logo_to_surface(surface);

  switch(kind) {
    case output_t::PDF:
      cairo_surface_flush(surface);
      break;
    case output_t::PNG:
      cairo_surface_write_to_png(surface, out_path.c_str());
      break;
  }
  cairo_surface_destroy (surface);
  return EXIT_SUCCESS;
}
};
