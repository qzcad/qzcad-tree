/**
  * @author Сергей Чопоров
  * @date 31/07/08
  * @version 1.0.0
  **/
#ifndef QSTDREDIRECTOR_H
#define QSTDREDIRECTOR_H

#include <iostream>
#include "qtxtsender.h"
/**
 * Класс-шаблон для перехвата стандартного потока сообщений
 * @see QTxtSender, std::streambuf
 */
template< class Elem = char, class Tr = std::char_traits< Elem > >
class QStdRedirector : public QTxtSender, public std::streambuf
{
public:
    /**
     * @brief Конструктор
     * @param a_Stream Ссылка на поток
     * @param parent Укзатель на родительский объект
     */
    QStdRedirector(std::ostream& a_Stream, QObject *parent = 0):
            QTxtSender(parent), m_Stream(a_Stream)
    {
        m_pBuf = m_Stream.rdbuf( this ); // sets new stream
    }
    /**
      * Деструктор
      */
    ~QStdRedirector()
    {
        m_Stream.rdbuf( m_pBuf ); // restores original stream
    }
    /**
     * @brief Реализация метода печати сообщения
     * @param _Ptr Указатель на сообщение
     * @param _Count Размер сообщения
     * @return Размер обработанного сообщения
     */
    std::streamsize xsputn( const Elem* _Ptr, std::streamsize _Count )
    {
        printMessage( _Ptr, _Count );
        return _Count;
    }
    /**
     * @brief Реализация обработки переполнения
     * @param v Информационный символ
     * @return Информационный символ
     */
    typename Tr::int_type overflow( typename Tr::int_type v )
    {
        Elem ch = Tr::to_char_type( v );
        printMessage( &ch, 1 );
        return Tr::not_eof( v );
    }

protected:
    std::basic_ostream<Elem, Tr>& m_Stream;
    std::streambuf*               m_pBuf;

};

#endif // QSTDREDIRECTOR_H
