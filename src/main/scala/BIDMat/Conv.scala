package BIDMat

object Conv {

  // 2d linear discrete convolution of a single-channel image with some filter.
  // ok maybe it's actually a 2d cross-correlation, but the neural network 
  // community is really bad at using proper terminology so whatever.
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

  def im2col_dim(image_width: Int, image_height: Int,
    filter_width: Int, filter_height: Int,
    pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int,
    dilation_width: Int, dilation_height: Int): (Int, Int, Int, Int) = {

      val out_width = (image_width + 2 * pad_width - (dilation_width * (filter_width - 1) + 1)) / stride_width + 1
      val out_height = (image_height + 2 * pad_height - (dilation_height * (filter_height - 1) + 1)) / stride_height + 1
      val out_cols = out_width * out_height
      val out_rows = filter_width * filter_height

      return (
          out_width,
          out_height,
          out_rows,
          out_cols
        )
  }

  def im2col(image: Array[Float], image_width: Int, image_height: Int, image_channels: Int,
    filter_width: Int, filter_height: Int, pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int, dilation_width: Int, dilation_height: Int,
    out: Array[Float]) = {

      val (out_width, out_height, out_rows, out_cols) = im2col_dim(image_width, image_height, filter_width, filter_height, pad_width, pad_height, stride_width, stride_height, dilation_width, dilation_height)
      val channel_size = image_width * image_height

      var channel_offset = 0
      var channel = 0
      var out_idx = 0
      while (channel < image_channels) {
        var out_col = 0
        var image_col_offset = -pad_width
        while (out_col < out_width) {
          var out_row = 0
          var image_row_offset = -pad_height
          while (out_row < out_height) {
            var filter_col_offset = 0
            while (filter_col_offset < filter_height) {
              var filter_row_offset = 0
              while (filter_row_offset < filter_width) {

                val image_row = image_row_offset + filter_row_offset
                val image_col = image_col_offset + filter_col_offset
                if (image_row >= 0 && image_row < image_height &&
                  image_col >= 0 && image_col < image_width) {
                    out(out_idx) = image(image_col * image_height + image_row)
                } else {
                    out(out_idx) = 0
                }
                out_idx += 1

                filter_row_offset += 1
              }
              filter_col_offset += 1
            }
            image_row_offset += stride_height
            out_row += 1
          }
          image_col_offset += stride_width
          out_col += 1
        }
        channel_offset += image_width * image_height
        channel += 1
      }

  }

  def flip_filter(filter: Array[Float], filter_width: Int, filter_height: Int) = {
    val filter_size = filter_width * filter_height
    var filter_idx = 0
    while (filter_idx < filter_size / 2) {
      val tmp = filter(filter_idx)
      filter(filter_idx) = filter(filter_size - filter_idx - 1)
      filter(filter_size - filter_idx - 1) = tmp
      filter_idx += 1
    }
  }

  def convolution_2d_mm(
    image: Array[Float], image_width: Int, image_height: Int, image_channels: Int,
    filter: Array[Float], filter_width: Int, filter_height: Int,
    pad_width: Int, pad_height: Int,
    stride_width: Int, stride_height: Int,
    dilation_width: Int, dilation_height: Int,
    out: Array[Float]): FMat = {

      val (out_width, out_height, out_rows, out_cols) = im2col_dim(image_width, image_height, filter_width, filter_height, pad_width, pad_height, stride_width, stride_height, dilation_width, dilation_height)

      im2col(image, image_width, image_height, image_channels,
        filter_width, filter_height,
        pad_width, pad_height,
        stride_width, stride_height,
        dilation_width, dilation_height,
        out)

      flip_filter(filter, filter_width, filter_height)
      val filter_mat = new FMat(1, filter_width * filter_height, filter)
      val conv_mat = new FMat(out_rows, out_cols, out)

      val out_mat = filter_mat * conv_mat
      return out_mat.view(out_width, out_height)
  }

}
