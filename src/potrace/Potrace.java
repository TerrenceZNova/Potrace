package potrace;


import org.opencv.core.Mat;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;


public class Potrace {

    public static boolean[][] bitmapToBinary(Mat srcImage) {

        assert (srcImage != null && !srcImage.empty());

        Mat binaryImage = new Mat();
        if (srcImage.channels() == 3) {
            Mat grayImage = new Mat();
            Imgproc.cvtColor(srcImage, grayImage, Imgproc.COLOR_BGR2GRAY);
            Imgproc.threshold(grayImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        } else {
            Imgproc.threshold(srcImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        }

        boolean[][] result = new boolean[srcImage.rows() + 2][srcImage.cols() + 2];
        for (int row = 0; row < binaryImage.rows(); row++) {
            for (int col = 0; col < binaryImage.cols(); col++) {
                if ((int) binaryImage.get(row, col)[0] == 255) {
                    result[row + 1][col + 1] = true;
                } else {
                    result[row + 1][col + 1] = false;
                }
            }
        }

        // white borders
        for (int col = 0; col < result[0].length; col++) {
            result[0][col] = true;
            result[result.length - 1][col] = true;
        }
        for (int row = 1; row < result.length - 1; row++) {
            result[row][0] = true;
            result[row][result[0].length - 1] = true;
        }

        return result;
    }


    //Below is for test
    public static void main(String[] args) {

        String dllPath = "C:\\OpenCV\\opencv\\build\\java\\x64\\opencv_java320.dll";
        System.load(dllPath);

        testBitmapToBinary();
    }

    private static void testBitmapToBinary() {

        String filePath = "E:\\Java_Projects\\Potrace\\resources\\sourceEntireImages\\4a.png";
        Mat srcImage = Imgcodecs.imread(filePath);
        boolean[][] result = bitmapToBinary(srcImage);

        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[i].length; j++) {
                if (result[i][j]) {
                    System.out.print(0);
                } else {
                    System.out.print(1);
                }
            }
            System.out.println();
        }
    }
}
