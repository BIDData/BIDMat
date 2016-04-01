package BIDMat

import MatFunctions._

object Conv {

  // direct 2d convolution
  def convolution_2d(image: Array[Float], image_width: Int, image_height: Int,
    filter: Array[Float], filter_width: Int, filter_height: Int,
    out: Array[Float],
    stride: Int = 1, pad_zeros: Boolean = true,
    image_offset: Int = 0, filter_offset: Int = 0, out_offset: Int = 0): Array[Float] = {

      // Compute output image dimensions
      val out_width = if (pad_zeros) {
          if (filter_width % 2 == 0) {
            image_width / stride + 1
          } else {
            image_width / stride
          }
        } else {
          (image_width - filter_width) / stride + 1
        }

      val out_height = if (pad_zeros) {
          if (filter_height % 2 == 0) {
            image_height / stride + 1
          } else {
            image_height / stride
          }
        } else {
          (image_height - filter_height) / stride + 1
        }

      val half_filter_width = filter_width >> 1
      val half_filter_height = filter_height >> 1

      // output indices
      var out_row = 0
      while (out_row < out_width) {
        var out_col = 0
        while (out_col < out_height) {
          var sum = 0.0f

          // filter indices
          var filter_row_offset = -half_filter_width
          while (filter_row_offset <= half_filter_width) {

            val image_row = if (pad_zeros) {
                out_row * stride + filter_row_offset
              } else {
                out_row * stride + filter_row_offset + half_filter_width
              }

            if (image_row >= 0 && image_row < image_width) {
              val filter_row = filter_width - (half_filter_width + filter_row_offset) - 1

              var filter_col_offset = -half_filter_height
              while (filter_col_offset <= half_filter_height) {

                val image_col = if (pad_zeros) {
                    out_col * stride + filter_col_offset
                  } else {
                    out_col * stride + filter_col_offset + half_filter_height
                  }

                val filter_col = filter_height - (half_filter_height + filter_col_offset) - 1

                if (image_col >= 0 && image_col < image_height) {
                  val image_idx = image_row * image_width + image_col + image_offset
                  val filter_idx = filter_row * filter_width + filter_col + filter_offset

                  sum += image(image_idx) * filter(filter_idx)
                }

                filter_col_offset += 1
              }
            }

            filter_row_offset += 1
          }

          val out_idx = out_row * out_width + out_col + out_offset
          out(out_idx) += sum

          out_col += stride
        }

        out_row += stride
      }

      out
  }

  // calculate im2col dimensions
  def im2col_dim(image: FND, filters: FND,
    pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int
    ): (Int, Int, Int, Int) = {

      val image_height = image.dim(0)
      val image_width = image.dim(1)
      val image_depth = image.dim(2)

      val filter_height = filters.dim(0)
      val filter_width = filters.dim(1)

      val out_width = (image_width + 2 * pad_width - filter_width) / stride_width + 1
      val out_height = (image_height + 2 * pad_height - filter_height) / stride_height + 1
      val out_rows = filter_width * filter_height * image_depth
      val out_cols = out_width * out_height

      return (
          out_width,
          out_height,
          out_rows,
          out_cols
        )
  }

  // perform im2col preprocessing step to lower spatial convolution to matrix multiplication
  def im2col(image: FND, filters: FND, 
    pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int) = {

      val image_height = image.dim(0)
      val image_width = image.dim(1)
      val image_depth = image.dim(2)

      val filter_height = filters.dim(0)
      val filter_width = filters.dim(1)

      val (out_width, out_height, out_rows, out_cols) = im2col_dim(
        image, filters, pad_width, pad_height, stride_width, stride_height)
      val image_size = image_width * image_height

      val out_mat = FMat(out_rows, out_cols)
      var out = out_mat.data

      var out_idx = 0
      var out_col = 0
      var image_col_offset = -pad_width
      while (out_col < out_width) {
        var out_row = 0
        var image_row_offset = -pad_height
        while (out_row < out_height) {

          var image_depth_offset = 0
          var depth = 0
          while (depth < image_depth) {

            var filter_col_offset = 0
            while (filter_col_offset < filter_height) {
              var filter_row_offset = 0
              while (filter_row_offset < filter_width) {

                val image_row = image_row_offset + filter_row_offset
                val image_col = image_col_offset + filter_col_offset
                if (image_row >= 0 && image_row < image_height &&
                  image_col >= 0 && image_col < image_width) {
                    out(out_idx) = image(image_col * image_height + image_row + image_depth_offset)
                } else {
                    out(out_idx) = 0
                }
                out_idx += 1

                filter_row_offset += 1
              }
              filter_col_offset += 1
            }
            image_depth_offset += image_size
            depth += 1
          }

          image_row_offset += stride_height
          out_row += 1
        }
        image_col_offset += stride_width
        out_col += 1
      }

      out_mat
  }

  // flip filters before multiplication so that we perform convolution and not
  // cross-correlation
  def flip_filters(filters: FND) = {
    val filter_height = filters.dim(0)
    val filter_width = filters.dim(1)
    val filter_depth = filters.dim(2)
    val num_filters = filters.dim(3)

    val filter_size = filter_width * filter_height
    var filter_layer = 0
    while (filter_layer < num_filters * filter_depth) {
      val filter_offset = filter_size * filter_layer
      var filter_idx = 0
      while (filter_idx < filter_size / 2) {
        val tmp = filters(filter_offset + filter_idx)
        filters(filter_offset + filter_idx) = filters(filter_offset + filter_size - filter_idx - 1)
        filters(filter_offset + filter_size - filter_idx - 1) = tmp
        filter_idx += 1
      }
      filter_layer += 1
    }
  }

  def stack_filters(filters: FND): FMat = {
    val filter_height = filters.dim(0)
    val filter_width = filters.dim(1)
    val filter_depth = filters.dim(2)
    val num_filters = filters.dim(3)
    val filter_vol = filter_height * filter_width * filter_depth

    var filter_mat = FMat(num_filters, filter_vol)
    var filter = 0
    while (filter < num_filters) {
      filter_mat(filter, ?) = filters(?, ?, ?, filter).toFMatView(1, filter_vol)
      filter += 1
    }
    filter_mat
  }

  // spatial convolution lowered to matrix multiplication
  def convolution_mm(
    image: FND, filters: FND,
    pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int): FND = {

      val filter_height = filters.dim(0)
      val filter_width = filters.dim(1)
      val num_filters = filters.dim(3)
      val image_depth = image.dim(2)
      val (out_width, out_height, out_rows, out_cols) = im2col_dim(
        image, filters, pad_width, pad_height, stride_width, stride_height)

      flip_filters(filters)
      val filter_mat = stack_filters(filters)

      val conv_mat = im2col(image, filters,
        pad_width, pad_height,
        stride_width, stride_height)

      val out_mat = filter_mat * conv_mat

      val out = FND(out_height, out_width, num_filters)
      var filter = 0
      while (filter < num_filters) {
        out(?, ?, filter) = FND(out_mat(filter, ?)).reshapeView(out_height, out_width)
        filter += 1
      }

      out
  }

}
