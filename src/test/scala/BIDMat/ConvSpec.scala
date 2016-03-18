package BIDMat

import BIDMat.MatFunctions._

class ConvSpec extends BIDMatSpec {

  "convolution_2d" should "work with default settings" in {
    val image = 1f\2\3\4 on
                5f\6\7\8 on
                9f\10\11\12 on
                13f\14\15\16

    val filter = 1f\1\1 on
                 1f\1\0 on
                 1f\0\0

    val out = zeros(4, 4)

    Conv.convolution_2d(image.data, image.nrows, image.ncols,
      filter.data, filter.nrows, filter.ncols, out.data)

    val expected = 14f\23\28\19 on
                   32\46\52\31 on
                   52\70\76\43 on
                   37\40\43\16

    assert_approx_eq(expected, out)
  }

  it should "with pad_zeros=false" in {
    val image = 1f\2\3\4 on
                5f\6\7\8 on
                9f\10\11\12 on
                13f\14\15\16

    val filter = 1f\1\1 on
                 1f\1\0 on
                 1f\0\0

    val out = zeros(2, 2)

    Conv.convolution_2d(image.data, image.nrows, image.ncols,
      filter.data, filter.nrows, filter.ncols, out.data,
      pad_zeros=false)

    val expected = 46f\52 on
                   70\76

    assert_approx_eq(expected, out)
  }

  "convolution_2d_mm" should "work" in {
    val image = 1f\2\3\4 on
                5f\6\7\8 on
                9f\10\11\12 on
                13f\14\15\16

    val filter = 1f\1\1 on
                 1f\1\0 on
                 1f\0\0

    val image_width = image.ncols
    val image_height = image.nrows
    val image_channels = 1
    val pad_width = 1
    val pad_height = 1
    val dilation_width = 1
    val dilation_height = 1
    val filter_width = filter.ncols
    val filter_height = filter.nrows
    val stride_width = 1
    val stride_height = 1
    val (out_width, out_height, out_rows, out_cols) = Conv.im2col_dim(image_width, image_height, filter_width, filter_height, pad_width, pad_height, stride_width, stride_height, dilation_width, dilation_height)

    val out = zeros(out_rows, out_cols)

    val actual = Conv.convolution_2d_mm(
        image.data, image_width, image_height, image_channels,
        filter.data, filter_width, filter_height,
        pad_width, pad_height,
        stride_width, stride_height,
        dilation_width, dilation_height,
        out.data
      )

    val expected = 14f\23\28\19 on
                   32\46\52\31 on
                   52\70\76\43 on
                   37\40\43\16

    assert_approx_eq(actual, expected)
  }

  it should "work without padding" in {
    val image = 1f\2\3\4 on
                5f\6\7\8 on
                9f\10\11\12 on
                13f\14\15\16

    val filter = 1f\1\1 on
                 1f\1\0 on
                 1f\0\0

    val image_width = image.ncols
    val image_height = image.nrows
    val image_channels = 1
    val pad_width = 0
    val pad_height = 0
    val dilation_width = 1
    val dilation_height = 1
    val filter_width = filter.ncols
    val filter_height = filter.nrows
    val stride_width = 1
    val stride_height = 1
    val (out_width, out_height, out_rows, out_cols) = Conv.im2col_dim(image_width, image_height, filter_width, filter_height, pad_width, pad_height, stride_width, stride_height, dilation_width, dilation_height)

    val out = zeros(out_rows, out_cols)

    val actual = Conv.convolution_2d_mm(
        image.data, image_width, image_height, image_channels,
        filter.data, filter_width, filter_height,
        pad_width, pad_height,
        stride_width, stride_height,
        dilation_width, dilation_height,
        out.data
      )

    val expected = 46f\52 on
                   70\76

    assert_approx_eq(actual, expected)
  }

}
