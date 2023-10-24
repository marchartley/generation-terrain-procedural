#ifndef SKELETONIZE_H
#define SKELETONIZE_H

// trace_skeleton.cpp
// Trace skeletonization result into polylines
//
// Lingdong Huang 2020

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <string>
#include <climits>

//================================
// ENUMS
//================================
#define HORIZONTAL 1
#define VERTICAL 2

//================================
// PARAMS
//================================
#define CHUNK_SIZE 10           // the chunk size
#define SAVE_RECTS 0            // additionally save bounding rects of chunks (for visualization)
#define MAX_ITER 9999            // maximum number of iterations


struct skeleton_tracer_t {
  //================================
  // GLOBALS
  //================================
  typedef unsigned char uchar;
  uchar* im; // the image
  int W;     // width
  int H;     // height

  skeleton_tracer_t(){
    im = NULL;
    rects.head = NULL;
    rects.tail = NULL;
  }

  //================================
  // DATASTRUCTURES
  //================================

  typedef struct _point_t {
    int x;
    int y;
    struct _point_t * next;
  } point_t;

  typedef struct _polyline_t {
    point_t* head;
    point_t* tail;
    struct _polyline_t* prev;
    struct _polyline_t* next;
    int size;
  } polyline_t;


  typedef struct _rect_t {
    int x;
    int y;
    int w;
    int h;
    struct _rect_t* next;
  } rect_t;

  struct _rects_t{
    rect_t* head;
    rect_t* tail;
  } rects;

  //================================
  // DATASTRUCTURE IMPLEMENTATION
  //================================

  polyline_t* new_polyline();
  std::string print_polyline(polyline_t* q);
  std::string print_polylines(polyline_t* q);
  void destroy_polylines(polyline_t* q);

  void reverse_polyline(polyline_t* q);

  void cat_tail_polyline(polyline_t* q0, polyline_t* q1);

  void cat_head_polyline(polyline_t* q0, polyline_t* q1);

  void add_point_to_polyline(polyline_t* q, int x, int y);

  polyline_t *prepend_polyline(polyline_t* q0, polyline_t* q1);

  std::string print_rects();

  void destroy_rects();

  void add_rect(int x, int y, int w, int h);

  //================================
  // RASTER SKELETONIZATION
  //================================
  // Binary image thinning (skeletonization) in-place.
  // Implements Zhang-Suen algorithm.
  // http://agcggs680.pbworks.com/f/Zhan-Suen_algorithm.pdf
  bool thinning_zs_iteration(int iter);;
  void thinning_zs();

  //================================
  // MAIN ALGORITHM
  //================================

  // check if a region has any white pixel
  int not_empty(int x, int y, int w, int h);

  /**merge ith fragment of second chunk to first chunk
   * @param c0   fragments from  first  chunk
   * @param c1i  ith fragment of second chunk
   * @param sx   (x or y) coordinate of the seam
   * @param isv  is vertical, not horizontal?
   * @param mode 2-bit flag,
   *             MSB = is matching the left (not right) end of the fragment from first  chunk
   *             LSB = is matching the right (not left) end of the fragment from second chunk
   * @return     matching successful?
   */
  int merge_impl(polyline_t* c0, polyline_t* c1i, int sx, int isv, int mode);

  /**merge fragments from two chunks
   * @param c0   fragments from first  chunk
   * @param c1   fragments from second chunk
   * @param sx   (x or y) coordinate of the seam
   * @param dr   merge direction, HORIZONTAL or VERTICAL?
   */
  polyline_t* merge_frags(polyline_t* c0, polyline_t* c1, int sx, int dr);

  /**recursive bottom: turn chunk into polyline fragments;
   * look around on 4 edges of the chunk, and identify the "outgoing" pixels;
   * add segments connecting these pixels to center of chunk;
   * apply heuristics to adjust center of chunk
   *
   * @param x    left of   chunk
   * @param y    top of    chunk
   * @param w    width of  chunk
   * @param h    height of chunk
   * @return     the polyline fragments
   */
  polyline_t* chunk_to_frags(int x, int y, int w, int h);


  /**Trace skeleton from thinning result.
   * Algorithm:
   * 1. if chunk size is small enough, reach recursive bottom and turn it into segments
   * 2. attempt to split the chunk into 2 smaller chunks, either horizontall or vertically;
   *    find the best "seam" to carve along, and avoid possible degenerate cases
   * 3. recurse on each chunk, and merge their segments
   *
   * @param x       left of   chunk
   * @param y       top of    chunk
   * @param w       width of  chunk
   * @param h       height of chunk
   * @param iter    current iteration
   * @return        an array of polylines
  */
  polyline_t* trace_skeleton(int x, int y, int w, int h, int iter);


  //================================
  // GUI/IO
  //================================
  void print_bitmap();

  char* trace(char* img, int w, int h);

  void destroy();

};


#endif // SKELETONIZE_H
