#include<iostream>
#include<fstream>

void create_files(int example_number)
{
    for (int i = 1; i <= example_number; i++)
    {
        std:: string path = "C:/Users/cilli/VS programmes/loop_file_testing/run_" + std:: to_string(i) + ".txt";
        std:: string path_ds = "C:/Users/cilli/VS programmes/loop_file_testing/run_ds" + std:: to_string(i) + ".txt"; 
        std:: fstream current_file(path,std::ios::out);
        std:: fstream current_file_ds(path_ds,std::ios::out);
        current_file.close();
        current_file_ds.close();
    }
}

int main()
{
    int x[100] = {1,2,3,4,5};
    std:: cout << sizeof(x) / 4;
    return 0;
}

void create_files()
{
    const char *path1,*path2,*path3,*path4,*path5,*path6,*path7,*path8,*path9,*path10,
    *path11,*path12,*path13,*path14,*path15,*path16,*path17,*path18,*path19,*path20,
    *path_ds1,*path_ds2,*path_ds3,*path_ds4,*path_ds5,*path_ds6,*path_ds7,*path_ds8,*path_ds9,*path_ds10,
    *path_ds11,*path_ds12,*path_ds13,*path_ds14,*path_ds15,*path_ds16,*path_ds17,*path_ds18,*path_ds19,*path_ds20;

    path1 = "C:/Users/cilli/VS programmes/Run_files1.0/run_1.txt";
    path2 = "C:/Users/cilli/VS programmes/Run_files1.0/run_2.txt";
    path3 = "C:/Users/cilli/VS programmes/Run_files1.0/run_3.txt";
    path4 = "C:/Users/cilli/VS programmes/Run_files1.0/run_4.txt";
    path5 = "C:/Users/cilli/VS programmes/Run_files1.0/run_5.txt";
    path6 = "C:/Users/cilli/VS programmes/Run_files1.0/run_6.txt";
    path7 = "C:/Users/cilli/VS programmes/Run_files1.0/run_7.txt";
    path8 = "C:/Users/cilli/VS programmes/Run_files1.0/run_8.txt";
    path9 = "C:/Users/cilli/VS programmes/Run_files1.0/run_9.txt";
    path10 = "C:/Users/cilli/VS programmes/Run_files1.0/run_10.txt";
    path11 = "C:/Users/cilli/VS programmes/Run_files1.0/run_11.txt";
    path12 = "C:/Users/cilli/VS programmes/Run_files1.0/run_12.txt";
    path13 = "C:/Users/cilli/VS programmes/Run_files1.0/run_13.txt";
    path14 = "C:/Users/cilli/VS programmes/Run_files1.0/run_14.txt";
    path15 = "C:/Users/cilli/VS programmes/Run_files1.0/run_15.txt";
    path16 = "C:/Users/cilli/VS programmes/Run_files1.0/run_16.txt";
    path17 = "C:/Users/cilli/VS programmes/Run_files1.0/run_17.txt";
    path18 = "C:/Users/cilli/VS programmes/Run_files1.0/run_18.txt";
    path19 = "C:/Users/cilli/VS programmes/Run_files1.0/run_19.txt";
    path20 = "C:/Users/cilli/VS programmes/Run_files1.0/run_20.txt";
    path_ds1 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds1.txt";
    path_ds2 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds2.txt";
    path_ds3 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds3.txt";
    path_ds4 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds4.txt";
    path_ds5 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds5.txt";
    path_ds6 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds6.txt";
    path_ds7 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds7.txt";
    path_ds8 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds8.txt";
    path_ds9 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds9.txt";
    path_ds10 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds10.txt";
    path_ds11 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds11.txt";
    path_ds12 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds12.txt";
    path_ds13 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds13.txt";
    path_ds14 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds14.txt";
    path_ds15 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds15.txt";
    path_ds16 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds16.txt";
    path_ds17 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds17.txt";
    path_ds18 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds18.txt";
    path_ds19 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds19.txt";
    path_ds20 = "C:/Users/cilli/VS programmes/Run_files1.0/run_ds20.txt";

    std::fstream run_1(path1,std::ios::out);
    std::fstream run_2(path2,std::ios::out);
    std::fstream run_3(path3,std::ios::out);
    std::fstream run_4(path4,std::ios::out);
    std::fstream run_5(path5,std::ios::out);
    std::fstream run_6(path6,std::ios::out);
    std::fstream run_7(path7,std::ios::out);
    std::fstream run_8(path8,std::ios::out);
    std::fstream run_9(path9,std::ios::out);
    std::fstream run_10(path10,std::ios::out);
    std::fstream run_11(path11,std::ios::out);
    std::fstream run_12(path12,std::ios::out);
    std::fstream run_13(path13,std::ios::out);
    std::fstream run_14(path14,std::ios::out);
    std::fstream run_15(path15,std::ios::out);
    std::fstream run_16(path16,std::ios::out);
    std::fstream run_17(path17,std::ios::out);
    std::fstream run_18(path18,std::ios::out);
    std::fstream run_19(path19,std::ios::out);
    std::fstream run_20(path20,std::ios::out);
    std::fstream run_ds1(path_ds1,std::ios::out);
    std::fstream run_ds2(path_ds2,std::ios::out);
    std::fstream run_ds3(path_ds3,std::ios::out);
    std::fstream run_ds4(path_ds4,std::ios::out);
    std::fstream run_ds5(path_ds5,std::ios::out);
    std::fstream run_ds6(path_ds6,std::ios::out);
    std::fstream run_ds7(path_ds7,std::ios::out);
    std::fstream run_ds8(path_ds8,std::ios::out);
    std::fstream run_ds9(path_ds9,std::ios::out);
    std::fstream run_ds10(path_ds10,std::ios::out);
    std::fstream run_ds11(path_ds11,std::ios::out);
    std::fstream run_ds12(path_ds12,std::ios::out);
    std::fstream run_ds13(path_ds13,std::ios::out);
    std::fstream run_ds14(path_ds14,std::ios::out);
    std::fstream run_ds15(path_ds15,std::ios::out);
    std::fstream run_ds16(path_ds16,std::ios::out);
    std::fstream run_ds17(path_ds17,std::ios::out);
    std::fstream run_ds18(path_ds18,std::ios::out);
    std::fstream run_ds19(path_ds19,std::ios::out);
    std::fstream run_ds20(path_ds20,std::ios::out);

    run_1.close(),run_2.close(),run_3.close(),run_4.close(),run_5.close(),
    run_6.close(),run_7.close(),run_8.close(),run_9.close(),run_10.close(),
    run_11.close(),run_12.close(),run_13.close(),run_14.close(),run_15.close(),
    run_16.close(),run_17.close(),run_18.close(),run_19.close(),run_20.close();
    run_ds1.close(),run_ds2.close(),run_ds3.close(),run_ds4.close(),run_ds5.close(),
    run_ds6.close(),run_ds7.close(),run_ds8.close(),run_ds9.close(),run_ds10.close(),
    run_ds11.close(),run_ds12.close(),run_ds13.close(),run_ds14.close(),run_ds15.close(),
    run_ds16.close(),run_ds17.close(),run_ds18.close(),run_ds19.close(),run_ds20.close();
}