#ifndef VIENNAMESH_CORE_ALGORITHM_PIPELINE_HPP
#define VIENNAMESH_CORE_ALGORITHM_PIPELINE_HPP

#include "viennamesh/core/algorithm.hpp"
#include "viennamesh/core/algorithm_factory.hpp"

#include "pugixml/pugixml.hpp"

namespace viennamesh
{

  class algorithm_pipeline
  {
  public:

    typedef std::map<string, algorithm_handle> AlgorithmMapType;


    bool add_algorithm( pugi::xml_node const & algorithm_node )
    {
      pugi::xml_attribute algorithm_name_attribute = algorithm_node.attribute("name");
      if ( algorithm_name_attribute.empty() )
      {
        error(1) << "Algorithm has no name attribute" << std::endl;
        return false;
      }
      string algorithm_name = algorithm_name_attribute.as_string();

      pugi::xml_attribute algorithm_type_attribute = algorithm_node.attribute("type");
      if ( algorithm_type_attribute.empty() )
      {
        error(1) << "Algorithm \"" << algorithm_name << "\" has no type attribute" << std::endl;
        return false;
      }
      string algorithm_type = algorithm_type_attribute.as_string();

      algorithm_handle algorithm = viennamesh::algorithm_factory().create_from_name( algorithm_type );

      if (!algorithm)
      {
        error(1) << "Algorithm creation from factory failed" << std::endl;
        return false;
      }

      for (pugi::xml_node paramater_node = algorithm_node.child("parameter");
           paramater_node;
           paramater_node = paramater_node.next_sibling("parameter"))
      {
        pugi::xml_attribute parameter_name_attribute = paramater_node.attribute("name");
        if (parameter_name_attribute.empty())
        {
          error(1) << "Parameter has no name attribute" << std::endl;
          return false;
        }
        string parameter_name = parameter_name_attribute.as_string();

        pugi::xml_attribute parameter_type_attribute = paramater_node.attribute("type");
        if (parameter_type_attribute.empty())
        {
          error(1) << "Parameter \"" << parameter_name << "\" has no type attribute" << std::endl;
          return false;
        }
        string parameter_type = parameter_type_attribute.as_string();


        string parameter_value = paramater_node.text().as_string();
        if (parameter_value.empty())
        {
          error(1) << "Parameter \"" << parameter_name << "\" has no value" << std::endl;
          return false;
        }


        if (parameter_type == "string")
        {
          algorithm->set_input( parameter_name, parameter_value );
        }
        else if (parameter_type == "double")
        {
          algorithm->set_input( parameter_name, lexical_cast<double>(parameter_value) );
        }
        else if (parameter_type == "point")
        {
          algorithm->set_input( parameter_name, stringtools::vector_from_string<double>(parameter_value) );
        }
        else if (parameter_type == "dynamic")
        {
          string source_algorithm_name = parameter_value.substr( 0, parameter_value.find("/") );
          string source_parameter_name = parameter_value.substr( parameter_value.find("/")+1 );

          std::map<string, algorithm_handle>::iterator ait = algorithms.find( source_algorithm_name );
          if (ait == algorithms.end())
          {
            error(1) << "Dynamic parameter \"" << parameter_name << "\": algorithm \"" << source_algorithm_name << "\" was not found" << std::endl;
            return false;
          }

          algorithm->set_input( parameter_name, make_parameter_link(ait->second, source_parameter_name) );
        }
        else
        {
          error(1) << "Parameter \"" << parameter_name << "\": type \"" << parameter_type << "\" is not supported" << std::endl;
          return false;
        }

      }


      std::pair<AlgorithmMapType::iterator, bool> result =
        algorithms.insert( std::make_pair(algorithm_name, algorithm) );

      if (!result.second)
      {
        error(1) << "Duplicated algorithm with \"" << algorithm_name << "\" and type \"" << algorithm_type << std::endl;
        return false;
      }

      algorithm_order.push_back( result.first );

      return true;
    }

    bool from_xml( pugi::xml_node const & xml )
    {
      for (pugi::xml_node algorithm_node = xml.child("algorithm");
           algorithm_node;
           algorithm_node = algorithm_node.next_sibling("algorithm"))
      {
        if (!add_algorithm( algorithm_node ))
        {
          clear();
          return false;
        }
      }

      return true;
    }

    bool run()
    {
      for (std::vector<AlgorithmMapType::iterator>::iterator ait = algorithm_order.begin(); ait != algorithm_order.end(); ++ait)
      {
        if (!(*ait)->second->run())
          return false;
      }

      return true;
    }

    void clear()
    {
      algorithms.clear();
    }

  private:

    AlgorithmMapType algorithms;
    std::vector<AlgorithmMapType::iterator> algorithm_order;
  };


}

#endif
