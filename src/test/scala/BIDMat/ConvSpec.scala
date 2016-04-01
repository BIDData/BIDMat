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

    assert_approx_eq(expected.data, out.data)
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

    assert_approx_eq(expected.data, out.data)
  }

  "convolution_mm" should "work" in {
    val image_c1 = 1f\2\3\4 on
                   5f\6\7\8 on
                   9f\10\11\12 on
                   13f\14\15\16

    val image_c2 = 3f\9\1\2 on
                   4f\9\1\0 on
                   2f\9\1\1 on
                   1f\9\1\2

    val image = FND(4, 4, 2)

    image(?,?,0) = FND(image_c1)
    image(?,?,1) = FND(image_c2)

    val filter1_c1 = 1f\1\1 on
                    1f\1\0 on
                    1f\0\0

    val filter1_c2 = -1f\0\1 on
                    0f\3\0 on
                    1f\0\ -1

    val filter2_c1 = 0f\1\1 on
                    1f\0\0 on
                    1f\0\ -1

    val filter2_c2 = 3f\4\1 on
                    0f\3\0 on
                    1f\4\3

    val filters = FND(3, 3, 2, 2)

    filters(?,?,0,0) = FND(filter1_c1)
    filters(?,?,1,0) = FND(filter1_c2)
    filters(?,?,0,1) = FND(filter2_c1)
    filters(?,?,1,1) = FND(filter2_c2)

    val pad_width = 1
    val pad_height = 1
    val stride_width = 1
    val stride_height = 1

    val actual = Conv.convolution_mm(
        image, filters,
        pad_width, pad_height,
        stride_width, stride_height)

    val expected_f1 = 14f\53\40\26 on
                      44\72\56\31 on
                      58\94\77\46 on
                      49\66\38\21

    val expected_f2 = 59f\84\33\22 on
                      85\142\83\36 on
                      91\156\96\39 on
                      44\87\53\2

    val expected = FND(4, 4, 2)
    expected(?,?,0) = FND(expected_f1)
    expected(?,?,1) = FND(expected_f2)

    assert_approx_eq(actual.data, expected.data)
  }

}
