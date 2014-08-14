/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.0
  * */
#ifndef MECHANICALPARAMETERS_H
#define MECHANICALPARAMETERS_H
/**
 * @brief Класс для передачи основных механических параметров (констант): модуль Юнга, коэффициент Пуассона
 */
class MechanicalParameters
{
public:
    MechanicalParameters();
    /**
     * @brief Конструктор для константных параметров
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     */
    MechanicalParameters(double E, double nu);
    /**
     * @brief Модуль Юнга
     * @return Значение модуля Юнга
     */
    double E() const;
    /**
     * @brief Установить значение модуля Юнга
     * @param E Значение модуля Юнга
     */
    void setE(double E);
    /**
     * @brief Коэффициент Пуассона
     * @return Значение коэффициента Пуассона
     */
    double nu() const;
    /**
     * @brief Установить значение коэффициента Пуассона
     * @param nu Значение коэффициента Пуассона
     */
    void setNu(double nu);
    /**
     * @brief Установить значения механических констант
     * @param E Значение модуля Юнга
     * @param nu Значение коэффициента Пуассона
     */
    void set(double E, double nu);
    /**
     * @brief Установить значение модуля Юнга и коэффициента Пуассона на основе значений модуля сдвига и объемного модуля упругости
     * @param G Модуль сдвига
     * @param K Объемный модуль упругости
     */
    void setByGK(double G, double K);
private:
    double E_; //!< Модуль Юнга
    double nu_; //!< Коэффициент Пуасссона
};

#endif // MECHANICALPARAMETERS_H
