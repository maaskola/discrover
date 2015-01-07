#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cairo.h>
#include <cairo-pdf.h>
#include <boost/lexical_cast.hpp>
#include "logo.hpp"

namespace Logo {
using namespace std;

enum class output_t { PDF, PNG };

const string font_face = "sans";
// const string font_face = "serif";

struct dimensions_t {
  dimensions_t(double scale)
      : factor(scale / 100.0),
        node_width(0.75 * scale),
        node_height(1.0 * scale),
        axis_rel_ext(0.1),
        axis_tick_len(15.0 * factor),
        axis_horiz_space(100.0 * factor),
        axis_offset(5.0 * factor),
        axis_line_width(2.0 * factor),
        axis_font_size(30.0 * factor),
        plot_margin(5.0 * factor){};

  double factor;
  double node_width;
  double node_height;
  double axis_rel_ext;
  double axis_tick_len;
  double axis_horiz_space;
  double axis_offset;
  double axis_line_width;
  double axis_font_size;
  double plot_margin;
};

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

const palette_t default_colors
    = {{0.0, 0.75, 0.0}, {0.0, 0.0, 1.0}, {1, 0.6470588, 0.0}, {1.0, 0.0, 0.0}};

// B58900 DC322F 6C71C4 2AA198
const palette_t solarized_colors = {{0xB5 / 255.0, 0x89 / 255.0, 0x00 / 255.0},
                                    {0xDC / 255.0, 0x32 / 255.0, 0x2F / 255.0},
                                    {0x6C / 255.0, 0x71 / 255.0, 0xC4 / 255.0},
                                    {0x2A / 255.0, 0xA1 / 255.0, 0x98 / 255.0}};

// 0FAD00 FF0000 FFC600 0011A4
const palette_t tetrad_colors = {{0x0F / 255.0, 0xAD / 255.0, 0x00 / 255.0},
                                 {0xFF / 255.0, 0x00 / 255.0, 0x00 / 255.0},
                                 {0xFF / 255.0, 0xC6 / 255.0, 0x00 / 255.0},
                                 {0x00 / 255.0, 0x11 / 255.0, 0xA4 / 255.0}};

using coords_t = vector<coord_t>;

coords_t letter_a = {{0, 0},
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
                     {3.6, 4}};

coords_t letter_t
    = {{10, 10}, {10, 9}, {6, 9}, {6, 0}, {4, 0}, {4, 9}, {0, 9}, {0, 10}};

void draw_letter_a_or_t(cairo_t *cr, const coords_t &coords,
                        const coord_t &origin, double width, double height,
                        rgb_t color) {
  cairo_move_to(cr, origin.x + coords[0].x / 10 * width,
                origin.y - coords[0].y / 10 * height);
  for (auto coord : coords)
    cairo_line_to(cr, origin.x + coord.x / 10 * width,
                  origin.y - coord.y / 10 * height);
  cairo_set_line_width(cr, 1.0);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_EVEN_ODD);
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_fill(cr);
}

void draw_letter_a(cairo_t *cr, const coord_t &coord, double width,
                   double height, rgb_t color) {
  draw_letter_a_or_t(cr, letter_a, coord, width, height, color);
}

void draw_letter_c(cairo_t *cr, const coord_t &coord, double width,
                   double height, rgb_t color) {
  cairo_save(cr);
  cairo_rectangle(cr, coord.x, coord.y - height, 0.5 * width, height);
  cairo_rectangle(cr, coord.x, coord.y - 0.35 * height, width, 0.35 * height);
  cairo_rectangle(cr, coord.x, coord.y - height, width, 0.35 * height);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);
  cairo_clip(cr);
  cairo_save(cr);
  cairo_translate(cr, coord.x + 0.5 * width, coord.y - 0.5 * height);
  cairo_scale(cr, 1.0 * width / height, 1.0);
  cairo_arc(cr, 0, 0, 0.5 * height, 0, M_PI * 2);
  cairo_arc(cr, 0, 0, 0.5 * 0.8 * height, 0, M_PI * 2);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_EVEN_ODD);
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_fill(cr);
  cairo_restore(cr);
  cairo_restore(cr);
}

void draw_letter_g(cairo_t *cr, const coord_t &coord, double width,
                   double height, rgb_t color) {
  draw_letter_c(cr, coord, width, height, color);
  cairo_rectangle(cr, coord.x + 0.5 * width, coord.y - 0.45 * height,
                  0.5 * width, 0.1 * height);
  cairo_rectangle(cr, coord.x + 0.9 * width, coord.y - 0.45 * height,
                  0.1 * width, 0.45 * height);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_fill(cr);
}

void draw_letter_t(cairo_t *cr, const coord_t &coord, double width,
                   double height, rgb_t color) {
  draw_letter_a_or_t(cr, letter_t, coord, width, height, color);
}

void draw_letter_u(cairo_t *cr, const coord_t &coord, double width,
                   double height, rgb_t color) {
  cairo_save(cr);
  cairo_rectangle(cr, coord.x, coord.y - 0.5 * height, width, 0.5 * height);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);
  cairo_clip(cr);
  cairo_save(cr);
  cairo_translate(cr, coord.x + 0.5 * width, coord.y - 0.5 * height);
  cairo_scale(cr, 1.0 * width / height, 1.0);
  cairo_arc(cr, 0, 0, 0.5 * height, 0, M_PI * 2);
  cairo_arc(cr, 0, 0, 0.5 * 0.8 * height, 0, M_PI * 2);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_EVEN_ODD);
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_fill(cr);
  cairo_restore(cr);
  cairo_restore(cr);
  cairo_rectangle(cr, coord.x, coord.y - height, 0.1 * width, 0.5 * height);
  cairo_rectangle(cr, coord.x + 0.9 * width, coord.y - height, 0.1 * width,
                  0.5 * height);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);
  cairo_set_source_rgb(cr, color.r, color.g, color.b);
  cairo_fill(cr);
}

void draw_letter(cairo_t *cr, size_t letter, const coord_t &coord, double width,
                 double height, const Options &options, double eps = 1e-6) {
  palette_t palette;
  switch (options.palette) {
    case Palette::Default:
      palette = default_colors;
      break;
    case Palette::Solarized:
      palette = solarized_colors;
      break;
    case Palette::Tetrad:
      palette = tetrad_colors;
      break;
  }
  if (height > eps)
    switch (letter) {
      case 0:
        draw_letter_a(cr, coord, width, height, palette.a);
        break;
      case 1:
        draw_letter_c(cr, coord, width, height, palette.c);
        break;
      case 2:
        draw_letter_g(cr, coord, width, height, palette.g);
        break;
      case 3:
        switch (options.alphabet) {
          case Alphabet::DNA:
            draw_letter_t(cr, coord, width, height, palette.t);
            break;
          case Alphabet::RNA:
            draw_letter_u(cr, coord, width, height, palette.t);
            break;
          case Alphabet::Undefined:
            if (options.revcomp)
              draw_letter_t(cr, coord, width, height, palette.t);
            else
              draw_letter_u(cr, coord, width, height, palette.t);
            break;
        }
        break;
      default:
        cout << "Error: wrong nucleotide index specified in generating "
                "sequence logo." << endl;
        exit(-1);
    }
}

double information_content(const column_t &col) {
  double x = 0;
  for (auto &y : col)
    if (y > 0)
      x -= y * log(y) / log(2.0);
  return 2 - x;
}

void draw_logo_to_surface(cairo_surface_t *surface, const matrix_t &matrix,
                          const dimensions_t &dims, const Options &options) {
  cairo_t *cr = cairo_create(surface);

  double basecol = 0;
  double baseline = dims.node_height;
  if (options.axes) {
    baseline *= 1.0 + dims.axis_rel_ext;
    baseline += dims.plot_margin;
    basecol = dims.axis_horiz_space;
  }
  coord_t current = {basecol, baseline};
  for (auto &col : matrix) {
    double col_height = 1;
    if (options.type == Type::Sequence)
      col_height = information_content(col) / 2;

    vector<size_t> order = {0, 1, 2, 3};
    if (options.order == Order::Frequency)
      sort(begin(order), end(order),
           [&](size_t a, size_t b) { return col[a] < col[b]; });

    for (auto idx : order) {
      double current_height = col[idx] * dims.node_height * col_height;
      draw_letter(cr, idx, current, dims.node_width, current_height, options);
      current.y -= current_height;
    }

    current.x += dims.node_width;
    current.y = baseline;
  }

  const double axis_basecol = basecol - dims.axis_offset;
  if (options.axes) {
    cairo_move_to(
        cr, axis_basecol,
        dims.node_height * (1.0 + 2 * dims.axis_rel_ext) + dims.plot_margin);
    cairo_line_to(cr, axis_basecol, dims.plot_margin);

    const vector<double> ticks = {0, 0.5, 1.0};
    for (auto pos : ticks) {
      cairo_move_to(cr, axis_basecol, dims.node_height * dims.axis_rel_ext
                                      + dims.node_height * pos
                                      + dims.plot_margin);
      cairo_rel_line_to(cr, -dims.axis_tick_len, 0);
    }

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, dims.axis_line_width);
    cairo_stroke(cr);

    const double axis_annot_basecol
        = axis_basecol - 1.75 * dims.axis_tick_len
          - (options.type == Type::Sequence ? 0 : 5);
    for (auto pos : ticks) {
      cairo_text_extents_t te;
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_select_font_face(cr, font_face.c_str(), CAIRO_FONT_SLANT_NORMAL,
                             CAIRO_FONT_WEIGHT_NORMAL);
      cairo_set_font_size(cr, dims.axis_font_size);
      string label = boost::lexical_cast<string>(
          (1.0 - pos) * (options.type == Type::Sequence ? 2 : 1));
      if (label == "0.5")
        label = "½";
      cairo_text_extents(cr, label.c_str(), &te);
      cairo_move_to(cr, axis_annot_basecol - te.width / 2 - te.x_bearing,
                    dims.node_height * dims.axis_rel_ext
                    + dims.node_height * pos - te.height / 2 - te.y_bearing
                    + dims.plot_margin);
      cairo_show_text(cr, label.c_str());
    }

    const double axis_label_basecol = axis_annot_basecol
                                      - 1.2 * dims.axis_font_size;
    cairo_save(cr);
    cairo_text_extents_t te;
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(cr, font_face.c_str(), CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, dims.axis_font_size);
    const string label = boost::lexical_cast<string>(
        options.type == Type::Sequence ? "IC [bit]" : "Freq.");
    cairo_text_extents(cr, label.c_str(), &te);
    cairo_move_to(cr, axis_label_basecol - te.height / 2 - te.y_bearing,
                  dims.node_height * dims.axis_rel_ext + dims.node_height * 0.5
                  + te.width / 2 + te.x_bearing);
    cairo_rotate(cr, -0.5 * M_PI);
    cairo_show_text(cr, label.c_str());
    cairo_restore(cr);
  }

  cairo_show_page(cr);
  cairo_destroy(cr);
}

string ending(output_t kind) {
  switch (kind) {
    case output_t::PDF:
      return "pdf";
    case output_t::PNG:
      return "png";
    default:
      return "";
  };
}

string draw_logo_sub(const matrix_t &matrix, const string &path, output_t kind,
                     const Options &options) {
  string out_path = path + "." + ending(kind);
  cout << "Sequence logo in " << out_path << endl;

  dimensions_t dims(options.scale);
  double width = dims.node_width * matrix.size();
  double height = dims.node_height;
  if (options.axes) {
    width += dims.axis_horiz_space;
    height *= 1.0 + 2 * dims.axis_rel_ext;
    height += 2 * dims.plot_margin;
  }

  cairo_surface_t *surface;
  switch (kind) {
    case output_t::PDF:
      surface = cairo_pdf_surface_create(out_path.c_str(), width, height);
      break;
    case output_t::PNG:
      surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
      break;
    default:
      return "";
  }

  draw_logo_to_surface(surface, matrix, dims, options);

  switch (kind) {
    case output_t::PDF:
      cairo_surface_flush(surface);
      break;
    case output_t::PNG:
      cairo_surface_write_to_png(surface, out_path.c_str());
      break;
  }
  cairo_surface_destroy(surface);
  return out_path;
}

vector<string> draw_logo_rc(const matrix_t &matrix, const string &path,
                            const Options &options) {
  vector<string> paths;
  if (options.pdf_logo)
    paths.push_back(draw_logo_sub(matrix, path, output_t::PDF, options));
  if (options.png_logo)
    paths.push_back(draw_logo_sub(matrix, path, output_t::PNG, options));
  return paths;
}

vector<string> draw_logo(const matrix_t &matrix, const string &out_path,
                         const Options &options) {
  if (not options.revcomp)
    return draw_logo_rc(matrix, out_path, options);
  else {
    vector<string> paths = draw_logo_rc(matrix, out_path + ".forward", options);
    auto rc_matrix = matrix;
    reverse(begin(rc_matrix), end(rc_matrix));
    for (auto &col : rc_matrix)
      reverse(begin(col), end(col));
    for (auto path : draw_logo_rc(rc_matrix, out_path + ".revcomp", options))
      paths.push_back(path);
    return paths;
  }
}
};
